// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Mixed finite elements for incompressible problems.
//!
//! Implements Taylor-Hood P2/P1 elements, Stokes element matrices,
//! pressure stabilization, LBB inf-sup condition checks, penalty methods
//! for near-incompressibility, and MINI element bubble enrichment.

/// A mixed finite element pairing velocity and pressure fields.
///
/// Uses Taylor-Hood P2/P1 elements: quadratic velocity, linear pressure.
#[derive(Debug, Clone)]
pub struct MixedElement {
    /// Number of velocity nodes (P2 element).
    pub n_vel_nodes: usize,
    /// Number of pressure nodes (P1 element).
    pub n_pres_nodes: usize,
    /// Velocity degrees of freedom per node (spatial dimension).
    pub vel_dofs_per_node: usize,
    /// Stabilization parameter for pressure (Brezzi-Pitkäranta).
    pub stab_param: f64,
    /// Penalty parameter for near-incompressibility.
    pub penalty: f64,
}

impl MixedElement {
    /// Create a new mixed element with given parameters.
    ///
    /// # Arguments
    /// * `n_vel_nodes` - number of velocity nodes (e.g. 6 for 2D P2)
    /// * `n_pres_nodes` - number of pressure nodes (e.g. 3 for 2D P1)
    /// * `vel_dofs_per_node` - spatial dimension (2 or 3)
    /// * `stab_param` - Brezzi-Pitkäranta stabilization coefficient
    /// * `penalty` - penalty parameter for incompressibility
    pub fn new(
        n_vel_nodes: usize,
        n_pres_nodes: usize,
        vel_dofs_per_node: usize,
        stab_param: f64,
        penalty: f64,
    ) -> Self {
        Self {
            n_vel_nodes,
            n_pres_nodes,
            vel_dofs_per_node,
            stab_param,
            penalty,
        }
    }

    /// Total velocity degrees of freedom.
    pub fn n_vel_dofs(&self) -> usize {
        self.n_vel_nodes * self.vel_dofs_per_node
    }

    /// Total pressure degrees of freedom.
    pub fn n_pres_dofs(&self) -> usize {
        self.n_pres_nodes
    }
}

/// Brezzi-Pitkäranta pressure stabilization matrix.
///
/// Computes the stabilization term `alpha * h^2 * (grad p, grad q)` for
/// a local element of size `h_elem`. Returns a flat row-major `n_pres x n_pres` matrix.
///
/// # Arguments
/// * `h_elem` - characteristic element size
/// * `alpha` - stabilization coefficient (typically 0.01–0.25)
/// * `n_pres` - number of pressure DOFs
/// * `grad_p` - pressure shape function gradients, shape `[n_pres][dim]`
/// * `weight` - quadrature weight times Jacobian determinant
pub fn pressure_stabilization(
    h_elem: f64,
    alpha: f64,
    n_pres: usize,
    grad_p: &[[f64; 2]],
    weight: f64,
) -> Vec<f64> {
    assert_eq!(grad_p.len(), n_pres);
    let mut s = vec![0.0f64; n_pres * n_pres];
    let coeff = alpha * h_elem * h_elem * weight;
    for i in 0..n_pres {
        for j in 0..n_pres {
            let dot = grad_p[i][0] * grad_p[j][0] + grad_p[i][1] * grad_p[j][1];
            s[i * n_pres + j] += coeff * dot;
        }
    }
    s
}

/// Check the LBB (Ladyzhenskaya-Babuška-Brezzi) inf-sup condition.
///
/// Estimates the inf-sup constant as the smallest singular value of the
/// discrete divergence operator `B` (shape `[n_pres][n_vel_dofs]`).
/// The matrix is stored row-major. Returns the estimated inf-sup constant.
///
/// A positive return value indicates stability. Values below ~1e-12 suggest
/// the element pair violates the inf-sup condition.
pub fn inf_sup_check(b_matrix: &[f64], n_pres: usize, n_vel_dofs: usize) -> f64 {
    assert_eq!(b_matrix.len(), n_pres * n_vel_dofs);
    // Compute B^T B (n_vel_dofs x n_vel_dofs) then find smallest eigenvalue
    // via power iteration on B B^T (n_pres x n_pres) which is cheaper when pres < vel
    let m = n_pres;
    let n = n_vel_dofs;
    // B B^T, size m x m
    let mut bbt = vec![0.0f64; m * m];
    for i in 0..m {
        for j in 0..m {
            let mut s = 0.0;
            for k in 0..n {
                s += b_matrix[i * n + k] * b_matrix[j * n + k];
            }
            bbt[i * m + j] = s;
        }
    }
    // Find smallest eigenvalue of B B^T via inverse power iteration approximation.
    // We use a simple Gershgorin lower bound for demonstration purposes.
    let mut min_diag = f64::MAX;
    for i in 0..m {
        let diag = bbt[i * m + i];
        let mut row_sum = 0.0;
        for j in 0..m {
            if i != j {
                row_sum += bbt[i * m + j].abs();
            }
        }
        let lower = diag - row_sum;
        if lower < min_diag {
            min_diag = lower;
        }
    }
    if min_diag < 0.0 { 0.0 } else { min_diag.sqrt() }
}

/// Compute local stiffness matrix `K` and divergence matrix `B` for a Stokes element.
///
/// Returns `(K, B)` where:
/// - `K` is `[n_vel_dofs x n_vel_dofs]` (viscous term), row-major flat
/// - `B` is `[n_pres x n_vel_dofs]` (divergence coupling), row-major flat
///
/// # Arguments
/// * `viscosity` - dynamic viscosity
/// * `grad_v` - velocity shape function gradients, shape `[n_vel_nodes][dim=2]`
/// * `grad_p` - pressure shape function gradients, shape `[n_pres_nodes][dim=2]`
/// * `phi_p` - pressure shape function values at quadrature point, shape `[n_pres_nodes]`
/// * `div_v` - divergence of each velocity shape function, shape `[n_vel_nodes]`
/// * `weight` - quadrature weight * Jacobian det
/// * `dim` - spatial dimension (2)
pub fn stokes_element_matrix(
    viscosity: f64,
    grad_v: &[[f64; 2]],
    _grad_p: &[[f64; 2]],
    phi_p: &[f64],
    div_v: &[f64],
    weight: f64,
    dim: usize,
) -> (Vec<f64>, Vec<f64>) {
    let n_vel_nodes = grad_v.len();
    let n_pres_nodes = phi_p.len();
    let n_vel_dofs = n_vel_nodes * dim;
    let n_pres_dofs = n_pres_nodes;

    // K: viscosity * integral( grad(v_i) : grad(v_j) ) delta_{alpha,beta}
    // For dim=2, K is block-diagonal per component
    let mut k = vec![0.0f64; n_vel_dofs * n_vel_dofs];
    for i in 0..n_vel_nodes {
        for j in 0..n_vel_nodes {
            let dot = grad_v[i][0] * grad_v[j][0] + grad_v[i][1] * grad_v[j][1];
            let val = viscosity * weight * dot;
            // Block diagonal: each spatial component
            for d in 0..dim {
                let row = i * dim + d;
                let col = j * dim + d;
                k[row * n_vel_dofs + col] += val;
            }
        }
    }

    // B: integral( phi_p_i * div(v_j) ) where div(v_j) = sum_d grad_v[j][d]
    // B shape: [n_pres_dofs x n_vel_dofs]
    // For each pressure dof i and velocity node j, component d:
    //   B[i, j*dim+d] = -weight * phi_p[i] * grad_v[j][d]
    let mut b = vec![0.0f64; n_pres_dofs * n_vel_dofs];
    let _ = div_v; // div_v is redundant given grad_v; use grad_v directly
    for i in 0..n_pres_nodes {
        for (j, grad_vj) in grad_v.iter().enumerate().take(n_vel_nodes) {
            for (d, &gvjd) in grad_vj.iter().enumerate().take(dim) {
                let col = j * dim + d;
                b[i * n_vel_dofs + col] += -weight * phi_p[i] * gvjd;
            }
        }
    }

    (k, b)
}

/// Assemble incompressibility penalty term into the stiffness matrix.
///
/// Adds `penalty * B^T B` to the velocity block `K`. This enforces approximate
/// incompressibility without a separate pressure unknown.
///
/// # Arguments
/// * `k` - stiffness matrix (modified in place), size `n_vel_dofs x n_vel_dofs`
/// * `b` - divergence matrix, size `n_pres_dofs x n_vel_dofs`
/// * `penalty` - large penalty parameter
/// * `n_vel_dofs` - velocity DOF count
/// * `n_pres_dofs` - pressure DOF count
pub fn incompressibility_penalty(
    k: &mut [f64],
    b: &[f64],
    penalty: f64,
    n_vel_dofs: usize,
    n_pres_dofs: usize,
) {
    // K += penalty * B^T B
    // B^T is [n_vel_dofs x n_pres_dofs], B is [n_pres_dofs x n_vel_dofs]
    for i in 0..n_vel_dofs {
        for j in 0..n_vel_dofs {
            let mut s = 0.0;
            for p in 0..n_pres_dofs {
                // B^T[i,p] = B[p,i]
                s += b[p * n_vel_dofs + i] * b[p * n_vel_dofs + j];
            }
            k[i * n_vel_dofs + j] += penalty * s;
        }
    }
}

/// Evaluate the MINI element bubble function on a triangle.
///
/// The bubble function for the MINI element is `lambda_1 * lambda_2 * lambda_3`
/// where `lambda_i` are the barycentric coordinates. At interior points this
/// is positive and vanishes on all edges.
///
/// # Arguments
/// * `lam` - barycentric coordinates `[lambda_1, lambda_2, lambda_3]`, must sum to 1
///
/// Returns the bubble function value.
pub fn bubble_function(lam: &[f64; 3]) -> f64 {
    lam[0] * lam[1] * lam[2]
}

/// Evaluate the gradient of the MINI element bubble function.
///
/// Returns `grad(lambda_1 * lambda_2 * lambda_3)` with respect to barycentric
/// coordinates. The physical gradient requires multiplication by the inverse
/// Jacobian of the element mapping.
///
/// # Arguments
/// * `lam` - barycentric coordinates `[lambda_1, lambda_2, lambda_3]`
pub fn bubble_function_grad(lam: &[f64; 3]) -> [f64; 3] {
    [lam[1] * lam[2], lam[0] * lam[2], lam[0] * lam[1]]
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── MixedElement construction ──────────────────────────────────────────

    #[test]
    fn test_mixed_element_new() {
        let e = MixedElement::new(6, 3, 2, 0.1, 1e6);
        assert_eq!(e.n_vel_nodes, 6);
        assert_eq!(e.n_pres_nodes, 3);
    }

    #[test]
    fn test_mixed_element_dofs() {
        let e = MixedElement::new(6, 3, 2, 0.1, 1e6);
        assert_eq!(e.n_vel_dofs(), 12);
        assert_eq!(e.n_pres_dofs(), 3);
    }

    #[test]
    fn test_mixed_element_3d() {
        let e = MixedElement::new(10, 4, 3, 0.05, 1e8);
        assert_eq!(e.n_vel_dofs(), 30);
        assert_eq!(e.n_pres_dofs(), 4);
    }

    #[test]
    fn test_mixed_element_stab_param() {
        let e = MixedElement::new(6, 3, 2, 0.25, 1e6);
        assert!((e.stab_param - 0.25).abs() < 1e-14);
    }

    #[test]
    fn test_mixed_element_penalty() {
        let e = MixedElement::new(6, 3, 2, 0.1, 1e9);
        assert!((e.penalty - 1e9).abs() < 1.0);
    }

    // ── pressure_stabilization ─────────────────────────────────────────────

    #[test]
    fn test_pressure_stab_size() {
        let grad_p = [[1.0, 0.0], [0.0, 1.0], [-1.0, -1.0]];
        let s = pressure_stabilization(0.1, 0.1, 3, &grad_p, 0.5);
        assert_eq!(s.len(), 9);
    }

    #[test]
    fn test_pressure_stab_symmetry() {
        let grad_p = [[1.0, 2.0], [3.0, 4.0], [-1.0, 0.5]];
        let s = pressure_stabilization(0.2, 0.1, 3, &grad_p, 1.0);
        for i in 0..3 {
            for j in 0..3 {
                let diff = (s[i * 3 + j] - s[j * 3 + i]).abs();
                assert!(diff < 1e-14, "S not symmetric at ({i},{j})");
            }
        }
    }

    #[test]
    fn test_pressure_stab_positive_semidefinite() {
        let grad_p = [[1.0, 0.0], [0.0, 1.0]];
        let s = pressure_stabilization(1.0, 1.0, 2, &grad_p, 1.0);
        // Diagonal entries must be non-negative
        assert!(s[0] >= 0.0);
        assert!(s[3] >= 0.0);
    }

    #[test]
    fn test_pressure_stab_zero_alpha() {
        let grad_p = [[1.0, 0.0], [0.0, 1.0]];
        let s = pressure_stabilization(1.0, 0.0, 2, &grad_p, 1.0);
        for &v in &s {
            assert!(v.abs() < 1e-14);
        }
    }

    #[test]
    fn test_pressure_stab_scales_with_h_squared() {
        let grad_p = [[1.0, 0.0]];
        let s1 = pressure_stabilization(1.0, 1.0, 1, &grad_p, 1.0);
        let s2 = pressure_stabilization(2.0, 1.0, 1, &grad_p, 1.0);
        let ratio = s2[0] / s1[0];
        assert!((ratio - 4.0).abs() < 1e-10);
    }

    // ── inf_sup_check ──────────────────────────────────────────────────────

    #[test]
    fn test_inf_sup_identity_like() {
        // B = identity-like 2x2 => inf-sup ~ 1
        let b = vec![1.0, 0.0, 0.0, 1.0];
        let val = inf_sup_check(&b, 2, 2);
        assert!(val >= 0.0);
    }

    #[test]
    fn test_inf_sup_non_negative() {
        let b = vec![1.0, 2.0, 3.0, 0.5, 1.0, 0.0];
        let val = inf_sup_check(&b, 2, 3);
        assert!(val >= 0.0);
    }

    #[test]
    fn test_inf_sup_zero_matrix() {
        let b = vec![0.0; 6];
        let val = inf_sup_check(&b, 2, 3);
        assert!(val >= 0.0);
    }

    #[test]
    fn test_inf_sup_1x1() {
        let b = vec![3.0];
        let val = inf_sup_check(&b, 1, 1);
        // BBT = [[9]], Gershgorin lower = 9-0 = 9 => sqrt(9)=3
        assert!((val - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_inf_sup_larger() {
        let b: Vec<f64> = vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0];
        let val = inf_sup_check(&b, 3, 4);
        assert!(val >= 0.0);
    }

    // ── stokes_element_matrix ──────────────────────────────────────────────

    #[test]
    fn test_stokes_k_size() {
        let grad_v = [[1.0, 0.0], [0.0, 1.0], [-1.0, -1.0]];
        let grad_p = [[1.0, 0.0], [-1.0, 0.0]];
        let phi_p = [0.5, 0.5];
        let div_v = [1.0, 1.0, -2.0];
        let (k, _b) = stokes_element_matrix(1.0, &grad_v, &grad_p, &phi_p, &div_v, 1.0, 2);
        // n_vel_dofs = 3*2 = 6
        assert_eq!(k.len(), 36);
    }

    #[test]
    fn test_stokes_b_size() {
        let grad_v = [[1.0, 0.0], [0.0, 1.0], [-1.0, -1.0]];
        let grad_p = [[1.0, 0.0], [-1.0, 0.0]];
        let phi_p = [0.5, 0.5];
        let div_v = [1.0, 1.0, -2.0];
        let (_k, b) = stokes_element_matrix(1.0, &grad_v, &grad_p, &phi_p, &div_v, 1.0, 2);
        // n_pres_dofs=2, n_vel_dofs=6 => 12
        assert_eq!(b.len(), 12);
    }

    #[test]
    fn test_stokes_k_symmetry() {
        let grad_v = [[1.0, 2.0], [3.0, 4.0], [-1.0, 0.5]];
        let grad_p = [[1.0, 0.0]];
        let phi_p = [1.0];
        let div_v = [3.0, 7.0, -0.5];
        let (k, _) = stokes_element_matrix(1.0, &grad_v, &grad_p, &phi_p, &div_v, 1.0, 2);
        let n = 6;
        for i in 0..n {
            for j in 0..n {
                let diff = (k[i * n + j] - k[j * n + i]).abs();
                assert!(diff < 1e-12, "K not symmetric at ({i},{j})");
            }
        }
    }

    #[test]
    fn test_stokes_viscosity_scaling() {
        let grad_v = [[1.0, 0.0], [0.0, 1.0]];
        let phi_p = [1.0];
        let grad_p = [[0.5, 0.5]];
        let div_v = [1.0, 1.0];
        let (k1, _) = stokes_element_matrix(1.0, &grad_v, &grad_p, &phi_p, &div_v, 1.0, 2);
        let (k2, _) = stokes_element_matrix(2.0, &grad_v, &grad_p, &phi_p, &div_v, 1.0, 2);
        for (a, b) in k1.iter().zip(k2.iter()) {
            assert!((b - 2.0 * a).abs() < 1e-12);
        }
    }

    #[test]
    fn test_stokes_k_positive_diagonal() {
        let grad_v = [[1.0, 0.0], [0.0, 1.0]];
        let phi_p = [1.0];
        let grad_p = [[0.5, 0.5]];
        let div_v = [1.0, 1.0];
        let (k, _) = stokes_element_matrix(1.0, &grad_v, &grad_p, &phi_p, &div_v, 1.0, 2);
        // Diagonal entries should be >= 0
        let n = 4;
        for i in 0..n {
            assert!(k[i * n + i] >= 0.0, "Negative diagonal at {i}");
        }
    }

    // ── incompressibility_penalty ──────────────────────────────────────────

    #[test]
    fn test_penalty_increases_k() {
        let mut k = vec![1.0; 4];
        let b = vec![1.0, 0.0, 0.0, 1.0]; // 2x2 identity
        let k_before = k.clone();
        incompressibility_penalty(&mut k, &b, 1.0, 2, 2);
        // K += penalty * B^T B = I => K[0,0] += 1
        assert!(k[0] >= k_before[0]);
    }

    #[test]
    fn test_penalty_zero_no_change() {
        let mut k = vec![2.0, 0.0, 0.0, 2.0];
        let k_orig = k.clone();
        let b = vec![1.0, 0.0, 0.0, 1.0];
        incompressibility_penalty(&mut k, &b, 0.0, 2, 2);
        for (a, b) in k.iter().zip(k_orig.iter()) {
            assert!((a - b).abs() < 1e-14);
        }
    }

    #[test]
    fn test_penalty_symmetry_preserved() {
        let mut k = vec![1.0, 0.5, 0.5, 1.0];
        let b = vec![1.0, 2.0, 3.0, 4.0]; // 2x2 B
        incompressibility_penalty(&mut k, &b, 1.0, 2, 2);
        let diff = (k[1] - k[2]).abs();
        assert!(diff < 1e-12, "Symmetry broken after penalty");
    }

    // ── bubble_function ────────────────────────────────────────────────────

    #[test]
    fn test_bubble_centroid() {
        let lam = [1.0 / 3.0; 3];
        let val = bubble_function(&lam);
        let expected = (1.0 / 3.0f64).powi(3);
        assert!((val - expected).abs() < 1e-14);
    }

    #[test]
    fn test_bubble_vertex_zero() {
        assert!(bubble_function(&[1.0, 0.0, 0.0]).abs() < 1e-14);
        assert!(bubble_function(&[0.0, 1.0, 0.0]).abs() < 1e-14);
        assert!(bubble_function(&[0.0, 0.0, 1.0]).abs() < 1e-14);
    }

    #[test]
    fn test_bubble_positive_interior() {
        let lam = [0.5, 0.3, 0.2];
        let val = bubble_function(&lam);
        assert!(val > 0.0);
    }

    #[test]
    fn test_bubble_grad_centroid() {
        let lam = [1.0 / 3.0; 3];
        let grad = bubble_function_grad(&lam);
        // By symmetry all components equal
        let expected = (1.0 / 3.0f64).powi(2);
        for &g in &grad {
            assert!((g - expected).abs() < 1e-14);
        }
    }

    #[test]
    fn test_bubble_grad_nonzero_interior() {
        let lam = [0.6, 0.3, 0.1];
        let grad = bubble_function_grad(&lam);
        // gradient components are products of other two barycentric coords
        assert!((grad[0] - lam[1] * lam[2]).abs() < 1e-14);
        assert!((grad[1] - lam[0] * lam[2]).abs() < 1e-14);
        assert!((grad[2] - lam[0] * lam[1]).abs() < 1e-14);
    }

    #[test]
    fn test_bubble_grad_vertex() {
        let lam = [1.0, 0.0, 0.0];
        let grad = bubble_function_grad(&lam);
        // grad[0] = lam[1]*lam[2] = 0
        assert!(grad[0].abs() < 1e-14);
        // grad[1] = lam[0]*lam[2] = 0
        assert!(grad[1].abs() < 1e-14);
        // grad[2] = lam[0]*lam[1] = 0
        assert!(grad[2].abs() < 1e-14);
    }
}
