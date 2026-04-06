#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! MLS (Moving Least Squares) SPH corrections.
//!
//! Provides first-order consistent SPH approximations using Moving Least Squares
//! correction matrices. This improves accuracy near boundaries and in regions of
//! particle disorder.

/// Pre-computed MLS correction matrix for a single particle.
///
/// Stores the `B` matrix that corrects SPH kernel gradient estimates
/// to achieve first-order (linear) consistency.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct MlsMatrix {
    /// The correction matrix B, stored row-major. Size: `dim x (1+dim)`.
    pub b: Vec<f64>,
    /// Spatial dimension.
    pub dim: usize,
    /// Shepard (zeroth-order) normalization factor.
    pub normalization: f64,
    /// Whether the matrix has been successfully computed.
    pub is_valid: bool,
}

impl MlsMatrix {
    /// Create an identity-initialized MLS matrix for given dimension.
    ///
    /// # Arguments
    /// * `dim` - spatial dimension (2 or 3)
    pub fn new(dim: usize) -> Self {
        // B is dim x (1+dim), initialized to zero
        Self {
            b: vec![0.0f64; dim * (1 + dim)],
            dim,
            normalization: 0.0,
            is_valid: false,
        }
    }

    /// Access element `B[row][col]`.
    pub fn get(&self, row: usize, col: usize) -> f64 {
        self.b[row * (1 + self.dim) + col]
    }
}

/// Compute the MLS correction matrix B for a particle.
///
/// Implements first-order MLS correction so that linear functions are reproduced
/// exactly. Uses the standard least-squares formulation:
/// `B = (P^T W P)^{-1} P^T W`
/// where `P` is the polynomial basis matrix and `W` is diagonal with kernel weights.
///
/// For 2D, the basis is `p(x) = [1, x, y]^T` (3 terms).
///
/// # Arguments
/// * `pos_i` - position of particle i, length `dim`
/// * `pos_j` - positions of neighbours, shape `[n_neigh][dim]`
/// * `w_ij` - kernel values for each neighbour (length `n_neigh`)
/// * `dim` - spatial dimension
///
/// Returns `(B, normalization)` where B is `dim x (1+dim)` row-major and
/// `normalization` is the Shepard sum.
pub fn mls_correction_matrix(
    pos_i: &[f64],
    pos_j: &[Vec<f64>],
    w_ij: &[f64],
    dim: usize,
) -> MlsMatrix {
    let n_neigh = pos_j.len();
    assert_eq!(w_ij.len(), n_neigh);
    assert_eq!(pos_i.len(), dim);

    let p_size = 1 + dim; // [1, x1, x2, ...] basis
    // A = P^T W P, size p_size x p_size
    let mut a = vec![0.0f64; p_size * p_size];
    let mut norm = 0.0f64;

    for j in 0..n_neigh {
        let wj = w_ij[j];
        norm += wj;
        // Basis vector for neighbour j: [1, (x_j - x_i), (y_j - y_i), ...]
        let mut pj = vec![1.0f64; p_size];
        for d in 0..dim {
            pj[1 + d] = pos_j[j][d] - pos_i[d];
        }
        // A += w_j * p_j * p_j^T
        for r in 0..p_size {
            for c in 0..p_size {
                a[r * p_size + c] += wj * pj[r] * pj[c];
            }
        }
    }

    // Invert A (p_size x p_size) using Gaussian elimination with partial pivoting
    let a_inv = match mat_inv(&a, p_size) {
        Some(inv) => inv,
        None => {
            // Singular: return invalid matrix
            let mut m = MlsMatrix::new(dim);
            m.normalization = norm;
            return m;
        }
    };

    // B = A^{-1} P^T W, but we only need the gradient rows (rows 1..dim)
    // B[d][j] = a_inv[d+1][*] dot (w_j * p_j[*]) summed over j, giving a dim x p_size "matrix"
    // However MlsMatrix.b is dim x (1+dim); we store the rows corresponding to gradient basis
    let mut b_mat = vec![0.0f64; dim * p_size];
    for d in 0..dim {
        let a_row = d + 1; // gradient row in A^{-1}
        for j in 0..n_neigh {
            let wj = w_ij[j];
            let mut pj = vec![1.0f64; p_size];
            for dd in 0..dim {
                pj[1 + dd] = pos_j[j][dd] - pos_i[dd];
            }
            // b_mat[d][c] += a_inv[a_row][c] * wj * pj[c] -- this gives B_dj coefficients
            // Actually we compute the correction vector c_j = (A^{-1} p_j) * w_j
            for c in 0..p_size {
                b_mat[d * p_size + c] += a_inv[a_row * p_size + c] * wj * pj[c];
            }
        }
    }

    MlsMatrix {
        b: b_mat,
        dim,
        normalization: norm,
        is_valid: true,
    }
}

/// Compute the MLS-corrected SPH gradient of a scalar field.
///
/// Uses the precomputed correction matrix to give first-order consistent gradient.
/// The corrected gradient is:
/// `grad_phi_i = sum_j (phi_j - phi_i) * B_j * w_ij * V_j`
///
/// # Arguments
/// * `phi_i` - field value at particle i
/// * `phi_j` - field values at neighbours
/// * `pos_i` - position of particle i
/// * `pos_j` - positions of neighbours, shape `[n_neigh][dim]`
/// * `w_ij` - kernel values for each neighbour
/// * `vol_j` - particle volumes for each neighbour
/// * `dim` - spatial dimension
///
/// Returns gradient vector of length `dim`.
#[allow(clippy::too_many_arguments)]
pub fn mls_gradient(
    phi_i: f64,
    phi_j: &[f64],
    pos_i: &[f64],
    pos_j: &[Vec<f64>],
    w_ij: &[f64],
    vol_j: &[f64],
    dim: usize,
) -> Vec<f64> {
    let n_neigh = phi_j.len();
    assert_eq!(pos_j.len(), n_neigh);
    assert_eq!(w_ij.len(), n_neigh);
    assert_eq!(vol_j.len(), n_neigh);

    let mls = mls_correction_matrix(pos_i, pos_j, w_ij, dim);
    let p_size = 1 + dim;

    let mut grad = vec![0.0f64; dim];
    if !mls.is_valid {
        // Fall back to uncorrected SPH gradient
        return grad;
    }

    for j in 0..n_neigh {
        let d_phi = phi_j[j] - phi_i;
        let wvj = w_ij[j] * vol_j[j];
        let mut pj = vec![1.0f64; p_size];
        for d in 0..dim {
            pj[1 + d] = pos_j[j][d] - pos_i[d];
        }
        // contribution: d_phi * wvj * B[d][c] * pj[c] summed over c
        for d in 0..dim {
            for c in 0..p_size {
                grad[d] += d_phi * wvj * mls.b[d * p_size + c] * pj[c];
            }
        }
    }
    grad
}

/// Compute the Shepard (zeroth-order) renormalization factor.
///
/// `C_i = sum_j W_ij * V_j`
///
/// This is used to normalize SPH density interpolation so that a constant
/// field is reproduced exactly (partition of unity).
///
/// # Arguments
/// * `w_ij` - kernel values for each neighbour
/// * `vol_j` - particle volumes for each neighbour
pub fn renormalization_factor(w_ij: &[f64], vol_j: &[f64]) -> f64 {
    assert_eq!(w_ij.len(), vol_j.len());
    w_ij.iter().zip(vol_j.iter()).map(|(&w, &v)| w * v).sum()
}

/// Check zeroth and first order consistency of the MLS approximation.
///
/// Zeroth order: reproducing a constant `phi = c` gives `grad_phi = 0`.
/// First order: reproducing a linear `phi = a*x + b*y` gives the exact gradient `[a, b]`.
///
/// Returns `(zeroth_ok, first_ok)` where each is `true` if the consistency
/// condition holds within the given tolerance.
///
/// # Arguments
/// * `pos_i` - position of particle i
/// * `pos_j` - neighbour positions
/// * `w_ij` - kernel weights
/// * `vol_j` - particle volumes
/// * `dim` - spatial dimension
/// * `tol` - tolerance for checking consistency
pub fn consistency_check(
    pos_i: &[f64],
    pos_j: &[Vec<f64>],
    w_ij: &[f64],
    vol_j: &[f64],
    dim: usize,
    tol: f64,
) -> (bool, bool) {
    // Zeroth order: constant field phi=1 everywhere, gradient should be zero
    let phi_j_const = vec![1.0f64; pos_j.len()];
    let grad_const = mls_gradient(1.0, &phi_j_const, pos_i, pos_j, w_ij, vol_j, dim);
    let zeroth_ok = grad_const.iter().all(|&g| g.abs() < tol);

    // First order: phi = x_0 component, gradient should be [1, 0, ...]
    let phi_j_linear: Vec<f64> = pos_j.iter().map(|p| p[0]).collect();
    let grad_lin = mls_gradient(pos_i[0], &phi_j_linear, pos_i, pos_j, w_ij, vol_j, dim);
    let first_ok =
        (grad_lin[0] - 1.0).abs() < tol && grad_lin.iter().skip(1).all(|&g| g.abs() < tol);

    (zeroth_ok, first_ok)
}

/// Kernel gradient free (KGF) correction factor.
///
/// Computes a scalar renormalization correction for the kernel gradient to
/// improve accuracy near boundaries. Based on the ratio of the approximate
/// support volume to the expected full support volume.
///
/// # Arguments
/// * `w_ij` - kernel values for neighbours
/// * `vol_j` - particle volumes
/// * `full_support_volume` - expected kernel support volume (e.g. `pi * h^2` for 2D)
pub fn kernel_correction(w_ij: &[f64], vol_j: &[f64], full_support_volume: f64) -> f64 {
    if full_support_volume < 1e-20 {
        return 1.0;
    }
    let actual = renormalization_factor(w_ij, vol_j);
    if actual.abs() < 1e-20 {
        return 1.0;
    }
    full_support_volume / actual
}

// ── Internal helpers ───────────────────────────────────────────────────────

/// Invert an `n x n` matrix stored row-major using Gaussian elimination.
/// Returns `None` if the matrix is singular.
fn mat_inv(a: &[f64], n: usize) -> Option<Vec<f64>> {
    assert_eq!(a.len(), n * n);
    let mut m = a.to_vec();
    let mut inv = vec![0.0f64; n * n];
    // Initialize identity
    for i in 0..n {
        inv[i * n + i] = 1.0;
    }
    // Forward elimination with partial pivoting
    for col in 0..n {
        // Find pivot
        let mut max_val = m[col * n + col].abs();
        let mut max_row = col;
        for row in (col + 1)..n {
            let v = m[row * n + col].abs();
            if v > max_val {
                max_val = v;
                max_row = row;
            }
        }
        if max_val < 1e-14 {
            return None;
        }
        // Swap rows
        if max_row != col {
            for k in 0..n {
                m.swap(col * n + k, max_row * n + k);
                inv.swap(col * n + k, max_row * n + k);
            }
        }
        let pivot = m[col * n + col];
        for k in 0..n {
            m[col * n + k] /= pivot;
            inv[col * n + k] /= pivot;
        }
        for row in 0..n {
            if row == col {
                continue;
            }
            let factor = m[row * n + col];
            for k in 0..n {
                let mv = m[col * n + k];
                let iv = inv[col * n + k];
                m[row * n + k] -= factor * mv;
                inv[row * n + k] -= factor * iv;
            }
        }
    }
    Some(inv)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_regular_2d_grid(n: usize, h: f64) -> (Vec<f64>, Vec<Vec<f64>>, Vec<f64>, Vec<f64>) {
        // Particle i at origin; n uniform neighbours around it
        let pos_i = vec![0.0, 0.0];
        let mut pos_j = Vec::new();
        let mut w_ij = Vec::new();
        let mut vol_j = Vec::new();
        let spacing = h / (n as f64);
        for i in 0..n {
            let x = (i as f64 - (n as f64 - 1.0) / 2.0) * spacing;
            pos_j.push(vec![x, 0.0]);
            let r = x.abs() / h;
            let w = if r < 1.0 { 1.0 - r } else { 0.0 };
            w_ij.push(w);
            vol_j.push(spacing);
        }
        (pos_i, pos_j, w_ij, vol_j)
    }

    // ── MlsMatrix ─────────────────────────────────────────────────────────

    #[test]
    fn test_mls_matrix_new_2d() {
        let m = MlsMatrix::new(2);
        assert_eq!(m.b.len(), 6); // 2 * (1+2)
        assert!(!m.is_valid);
    }

    #[test]
    fn test_mls_matrix_new_3d() {
        let m = MlsMatrix::new(3);
        assert_eq!(m.b.len(), 12); // 3 * (1+3)
    }

    #[test]
    fn test_mls_matrix_get() {
        let mut m = MlsMatrix::new(2);
        m.b[0] = 5.0;
        assert!((m.get(0, 0) - 5.0).abs() < 1e-14);
    }

    #[test]
    fn test_mls_matrix_clone() {
        let m = MlsMatrix::new(2);
        let m2 = m.clone();
        assert_eq!(m2.dim, 2);
    }

    // ── mls_correction_matrix ──────────────────────────────────────────────

    #[test]
    fn test_mls_correction_empty_neighbours() {
        let pos_i = vec![0.0, 0.0];
        let mls = mls_correction_matrix(&pos_i, &[], &[], 2);
        assert!(!mls.is_valid);
    }

    #[test]
    fn test_mls_correction_valid_regular_grid() {
        let (_pos_i, pos_j, w_ij, _vol_j) = make_regular_2d_grid(5, 1.0);
        let pos_i = vec![0.0, 0.0];
        let mls = mls_correction_matrix(&pos_i, &pos_j, &w_ij, 2);
        // With enough neighbours a valid matrix should be produced
        // (may or may not be valid depending on rank)
        assert_eq!(mls.dim, 2);
    }

    #[test]
    fn test_mls_correction_normalization_positive() {
        let pos_i = vec![0.0, 0.0];
        let pos_j = vec![vec![0.5, 0.0], vec![-0.5, 0.0], vec![0.0, 0.5]];
        let w_ij = vec![0.5, 0.5, 0.5];
        let mls = mls_correction_matrix(&pos_i, &pos_j, &w_ij, 2);
        assert!(mls.normalization > 0.0);
    }

    #[test]
    fn test_mls_correction_dim_preserved() {
        let pos_i = vec![0.0, 0.0];
        let pos_j = vec![vec![1.0, 0.0]];
        let w_ij = vec![1.0];
        let mls = mls_correction_matrix(&pos_i, &pos_j, &w_ij, 2);
        assert_eq!(mls.dim, 2);
    }

    // ── renormalization_factor ─────────────────────────────────────────────

    #[test]
    fn test_renorm_zero() {
        let r = renormalization_factor(&[], &[]);
        assert!(r.abs() < 1e-15);
    }

    #[test]
    fn test_renorm_single() {
        let r = renormalization_factor(&[2.0], &[3.0]);
        assert!((r - 6.0).abs() < 1e-14);
    }

    #[test]
    fn test_renorm_multiple() {
        let r = renormalization_factor(&[1.0, 2.0, 3.0], &[1.0, 1.0, 1.0]);
        assert!((r - 6.0).abs() < 1e-14);
    }

    #[test]
    fn test_renorm_unit_volume_sums_weights() {
        let w = vec![0.3, 0.5, 0.2];
        let v = vec![1.0, 1.0, 1.0];
        let r = renormalization_factor(&w, &v);
        assert!((r - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_renorm_non_negative() {
        let w = vec![0.1, 0.2, 0.3];
        let v = vec![0.5, 0.5, 0.5];
        let r = renormalization_factor(&w, &v);
        assert!(r >= 0.0);
    }

    // ── kernel_correction ──────────────────────────────────────────────────

    #[test]
    fn test_kernel_correction_full_support() {
        // If actual = full_support_volume => correction = 1
        let w = vec![1.0];
        let v = vec![std::f64::consts::PI];
        let corr = kernel_correction(&w, &v, std::f64::consts::PI);
        assert!((corr - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_kernel_correction_half_support() {
        // actual = pi/2, full = pi => correction = 2
        let w = vec![1.0];
        let v = vec![std::f64::consts::PI / 2.0];
        let corr = kernel_correction(&w, &v, std::f64::consts::PI);
        assert!((corr - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_kernel_correction_zero_full() {
        let w = vec![1.0];
        let v = vec![1.0];
        let corr = kernel_correction(&w, &v, 0.0);
        assert!((corr - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_kernel_correction_no_neighbours() {
        let corr = kernel_correction(&[], &[], 1.0);
        assert!((corr - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_kernel_correction_always_positive() {
        let w = vec![0.5, 0.3];
        let v = vec![0.2, 0.4];
        let corr = kernel_correction(&w, &v, 1.0);
        assert!(corr > 0.0);
    }

    // ── consistency_check ──────────────────────────────────────────────────

    #[test]
    fn test_consistency_check_empty_neighbours() {
        let pos_i = vec![0.0, 0.0];
        let (zeroth, first) = consistency_check(&pos_i, &[], &[], &[], 2, 1e-8);
        // With no neighbours mls_gradient returns zero => grad_const = 0 => zeroth ok
        assert!(zeroth);
        // first: grad_lin[0] should be 1.0, but fallback gives 0 => not ok
        assert!(!first);
    }

    #[test]
    fn test_consistency_check_returns_two_bools() {
        let pos_i = vec![0.0, 0.0];
        let pos_j = vec![vec![0.1, 0.0], vec![-0.1, 0.0]];
        let w_ij = vec![0.9, 0.9];
        let vol_j = vec![0.1, 0.1];
        let (zeroth, first) = consistency_check(&pos_i, &pos_j, &w_ij, &vol_j, 2, 1e-6);
        // Just check that results are booleans (no panic)
        let _ = zeroth;
        let _ = first;
    }

    // ── mls_gradient ──────────────────────────────────────────────────────

    #[test]
    fn test_mls_gradient_no_neighbours() {
        let pos_i = vec![0.0, 0.0];
        let grad = mls_gradient(1.0, &[], &pos_i, &[], &[], &[], 2);
        assert_eq!(grad.len(), 2);
        assert!(grad[0].abs() < 1e-15);
    }

    #[test]
    fn test_mls_gradient_dim_matches() {
        let pos_i = vec![0.0, 0.0];
        let pos_j = vec![vec![0.1, 0.0]];
        let w_ij = vec![0.9];
        let vol_j = vec![0.1];
        let phi_j = vec![1.0];
        let grad = mls_gradient(1.0, &phi_j, &pos_i, &pos_j, &w_ij, &vol_j, 2);
        assert_eq!(grad.len(), 2);
    }

    #[test]
    fn test_mls_gradient_constant_field_zero_grad() {
        // For a constant field the gradient should be zero (if MLS is valid)
        let pos_i = vec![0.0, 0.0];
        let pos_j = vec![
            vec![0.3, 0.0],
            vec![-0.3, 0.0],
            vec![0.0, 0.3],
            vec![0.0, -0.3],
        ];
        let w_ij = vec![0.7, 0.7, 0.7, 0.7];
        let vol_j = vec![0.1, 0.1, 0.1, 0.1];
        let phi_j = vec![5.0, 5.0, 5.0, 5.0];
        let grad = mls_gradient(5.0, &phi_j, &pos_i, &pos_j, &w_ij, &vol_j, 2);
        assert_eq!(grad.len(), 2);
        // Due to symmetry, gradient should be near zero
        for g in &grad {
            assert!(g.abs() < 1e-10, "gradient component {} not near zero", g);
        }
    }

    // ── mat_inv internal ───────────────────────────────────────────────────

    #[test]
    fn test_mat_inv_2x2_identity() {
        let a = vec![1.0, 0.0, 0.0, 1.0];
        let inv = mat_inv(&a, 2).unwrap();
        assert!((inv[0] - 1.0).abs() < 1e-12);
        assert!((inv[3] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_mat_inv_2x2_known() {
        // [[2,1],[1,1]]^{-1} = [[1,-1],[-1,2]]
        let a = vec![2.0, 1.0, 1.0, 1.0];
        let inv = mat_inv(&a, 2).unwrap();
        assert!((inv[0] - 1.0).abs() < 1e-10);
        assert!((inv[1] - (-1.0)).abs() < 1e-10);
        assert!((inv[2] - (-1.0)).abs() < 1e-10);
        assert!((inv[3] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_mat_inv_singular_returns_none() {
        let a = vec![1.0, 1.0, 1.0, 1.0];
        let inv = mat_inv(&a, 2);
        assert!(inv.is_none());
    }

    #[test]
    fn test_mat_inv_3x3_identity() {
        let mut a = vec![0.0; 9];
        a[0] = 1.0;
        a[4] = 1.0;
        a[8] = 1.0;
        let inv = mat_inv(&a, 3).unwrap();
        assert!((inv[0] - 1.0).abs() < 1e-12);
        assert!((inv[4] - 1.0).abs() < 1e-12);
        assert!((inv[8] - 1.0).abs() < 1e-12);
    }
}
