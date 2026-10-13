// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! MLS (Moving Least Squares) SPH corrections.
//!
//! Provides first-order consistent SPH approximations using Moving Least Squares
//! correction matrices. This improves accuracy near boundaries and in regions of
//! particle disorder.

/// Pre-computed MLS correction matrix for a single particle.
///
/// For a particle `i` with neighbours `j` the moving-least-squares fit of a
/// scalar field uses the polynomial basis `p_j = [1, x_j - x_i]` and the moment
/// matrix `A = Σ_j w_j p_j p_jᵀ`.  This type stores the *gradient rows* of the
/// inverse moment matrix, `B[d][c] = (A⁻¹)[d + 1][c]`, so that the
/// first-order consistent correction vector of neighbour `j` is
/// `(B p_j)[d]` and
///
/// `∂φ/∂x_d ≈ Σ_j w_j (B p_j)[d] (φ_j - φ_i)`
///
/// reproduces every linear field exactly.
#[derive(Debug, Clone)]
pub struct MlsMatrix {
    /// Gradient rows of `A⁻¹`, stored row-major. Size: `dim x (1+dim)`.
    pub b: Vec<f64>,
    /// Spatial dimension.
    pub dim: usize,
    /// Shepard (zeroth-order) normalization factor `Σ_j w_j`.
    pub normalization: f64,
    /// Whether the matrix has been successfully computed.
    pub is_valid: bool,
}

impl MlsMatrix {
    /// Create a zero-initialized, invalid MLS matrix for given dimension.
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
/// Builds the moment matrix `A = Σ_j w_j p_j p_jᵀ` of the weighted
/// least-squares fit with basis `p_j = [1, x_j - x_i]` (3 terms in 2-D, 4 in
/// 3-D), inverts it and stores the gradient rows `B[d][c] = (A⁻¹)[d + 1][c]`
/// (see [`MlsMatrix`]).
///
/// The weights are the *quadrature* weights of the neighbours.  For the
/// volume-weighted fit used by [`mls_gradient`] pass `w_ij * V_j`; if only the
/// kernel values are passed the correction is exact for the unit-volume
/// quadrature `V_j = 1`.
///
/// The moment matrix is non-dimensionalised with the largest neighbour
/// distance and the total weight before it is inverted, so the singularity
/// test does not depend on the length unit (a nano-metre stencil is treated
/// like a metre-sized one).
///
/// # Arguments
/// * `pos_i` - position of particle i, length `dim`
/// * `pos_j` - positions of neighbours, shape `[n_neigh][dim]`
/// * `w_ij` - quadrature weight for each neighbour (length `n_neigh`)
/// * `dim` - spatial dimension
///
/// Returns the matrix; `is_valid` is `false` when the inputs have
/// inconsistent lengths or the neighbours do not span the `dim`-dimensional
/// space (e.g. fewer than `dim + 1` non-degenerate neighbours), in which case
/// `b` is zero and only `normalization` is meaningful.
pub fn mls_correction_matrix(
    pos_i: &[f64],
    pos_j: &[Vec<f64>],
    w_ij: &[f64],
    dim: usize,
) -> MlsMatrix {
    let mut result = MlsMatrix::new(dim);
    let n_neigh = pos_j.len();
    if w_ij.len() != n_neigh || pos_i.len() != dim || pos_j.iter().any(|p| p.len() != dim) {
        return result;
    }
    result.normalization = w_ij.iter().sum();

    // Characteristic length (largest distance of a contributing neighbour) and
    // total weight, used to non-dimensionalise the moment matrix.
    let mut length = 0.0f64;
    let mut weight_scale = 0.0f64;
    for (pj, &wj) in pos_j.iter().zip(w_ij.iter()) {
        if wj == 0.0 {
            continue;
        }
        let dist2: f64 = (0..dim).map(|d| (pj[d] - pos_i[d]).powi(2)).sum();
        length = length.max(dist2.sqrt());
        weight_scale += wj.abs();
    }
    if !(length > 0.0 && length.is_finite() && weight_scale > 0.0 && weight_scale.is_finite()) {
        return result;
    }

    let p_size = 1 + dim; // [1, x1, x2, ...] basis
    // Scaled moment matrix A_s = Σ_j (w_j / W) q_j q_jᵀ with q_j = [1, (x_j - x_i) / L].
    let mut a_scaled = vec![0.0f64; p_size * p_size];
    let mut q = vec![1.0f64; p_size];
    for (pj, &wj) in pos_j.iter().zip(w_ij.iter()) {
        let weight = wj / weight_scale;
        for d in 0..dim {
            q[1 + d] = (pj[d] - pos_i[d]) / length;
        }
        for r in 0..p_size {
            for c in 0..p_size {
                a_scaled[r * p_size + c] += weight * q[r] * q[c];
            }
        }
    }

    // Invert A_s with Gaussian elimination and partial pivoting.
    let Some(inv_scaled) = mat_inv(&a_scaled, p_size) else {
        return result;
    };

    // A = W S A_s S with S = diag(1, L, ..., L)  =>  A⁻¹ = S⁻¹ A_s⁻¹ S⁻¹ / W.
    // Keep only the gradient rows (d + 1) of A⁻¹.
    for d in 0..dim {
        let row = d + 1;
        for c in 0..p_size {
            let column_scale = if c == 0 { 1.0 } else { length };
            result.b[d * p_size + c] =
                inv_scaled[row * p_size + c] / (weight_scale * length * column_scale);
        }
    }
    if result.b.iter().any(|v| !v.is_finite()) {
        result.b.iter_mut().for_each(|v| *v = 0.0);
        return result;
    }
    result.is_valid = true;
    result
}

/// Compute the MLS-corrected SPH gradient of a scalar field.
///
/// With the volume-weighted quadrature weights `wv_j = w_ij V_j` and the
/// moment matrix `A = Σ_j wv_j p_j p_jᵀ` (`p_j = [1, x_j - x_i]`), the
/// first-order consistent gradient is
///
/// `∂φ/∂x_d = Σ_j wv_j (A⁻¹ p_j)[d + 1] (φ_j - φ_i)`
///
/// which is exact for every linear field `φ` because `A⁻¹ A = I`.  A constant
/// field gives a zero gradient by construction.
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
/// Returns the gradient vector of length `dim`, or `None` when the gradient is
/// not determined: the slices have inconsistent lengths or the neighbours do
/// not span the `dim`-dimensional space (fewer than `dim + 1` non-degenerate
/// neighbours, e.g. at an isolated particle).  No kernel gradient is available
/// here, so there is no uncorrected SPH gradient to fall back to; the caller
/// decides how to handle such particles.
pub fn mls_gradient(
    phi_i: f64,
    phi_j: &[f64],
    pos_i: &[f64],
    pos_j: &[Vec<f64>],
    w_ij: &[f64],
    vol_j: &[f64],
    dim: usize,
) -> Option<Vec<f64>> {
    let n_neigh = phi_j.len();
    if pos_j.len() != n_neigh || w_ij.len() != n_neigh || vol_j.len() != n_neigh {
        return None;
    }

    let weights: Vec<f64> = w_ij.iter().zip(vol_j.iter()).map(|(w, v)| w * v).collect();
    let mls = mls_correction_matrix(pos_i, pos_j, &weights, dim);
    if !mls.is_valid {
        return None;
    }
    let p_size = 1 + dim;

    let mut grad = vec![0.0f64; dim];
    for j in 0..n_neigh {
        let d_phi = phi_j[j] - phi_i;
        let wvj = weights[j];
        // (A⁻¹ p_j)[d + 1] = B[d][0] + Σ_m B[d][1 + m] (x_j - x_i)[m]
        for (d, gd) in grad.iter_mut().enumerate() {
            let mut correction = mls.b[d * p_size];
            for m in 0..dim {
                correction += mls.b[d * p_size + 1 + m] * (pos_j[j][m] - pos_i[m]);
            }
            *gd += wvj * correction * d_phi;
        }
    }
    if grad.iter().any(|g| !g.is_finite()) {
        return None;
    }
    Some(grad)
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
/// First order: reproducing a linear `phi = x_d` gives the exact gradient `e_d`
/// for every coordinate direction `d` (by linearity this covers every linear
/// field `a·x + b`).
///
/// Returns `(zeroth_ok, first_ok)` where each is `true` if the consistency
/// condition holds within the given tolerance.  If the MLS system is singular
/// (see [`mls_gradient`]) no consistency can be established and both flags are
/// `false`.
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
    let Some(grad_const) = mls_gradient(1.0, &phi_j_const, pos_i, pos_j, w_ij, vol_j, dim) else {
        return (false, false);
    };
    let zeroth_ok = grad_const.iter().all(|&g| g.abs() < tol);

    // First order: phi = x_d for each axis, gradient should be the unit vector e_d
    let mut first_ok = true;
    for axis in 0..dim {
        let phi_j_linear: Vec<f64> = pos_j.iter().map(|p| p[axis]).collect();
        let Some(grad_lin) =
            mls_gradient(pos_i[axis], &phi_j_linear, pos_i, pos_j, w_ij, vol_j, dim)
        else {
            return (false, false);
        };
        first_ok &= grad_lin.iter().enumerate().all(|(d, &g)| {
            let expected = if d == axis { 1.0 } else { 0.0 };
            (g - expected).abs() < tol
        });
    }

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

    /// Disordered 2-D stencil: particle i away from the origin with 8
    /// irregularly placed neighbours, non-uniform kernel values and
    /// non-uniform volumes.  Returns `(pos_i, pos_j, w_ij, vol_j)`.
    fn disordered_stencil_2d() -> (Vec<f64>, Vec<Vec<f64>>, Vec<f64>, Vec<f64>) {
        let pos_i = vec![0.3, -0.2];
        let offsets = [
            [0.11, 0.02],
            [-0.09, 0.05],
            [0.03, -0.12],
            [-0.04, 0.10],
            [0.08, 0.09],
            [-0.12, -0.07],
            [0.06, -0.05],
            [-0.02, -0.03],
        ];
        let pos_j = offsets
            .iter()
            .map(|o| vec![pos_i[0] + o[0], pos_i[1] + o[1]])
            .collect();
        let w_ij = vec![0.9, 0.7, 0.6, 0.8, 0.5, 0.4, 1.1, 1.3];
        let vol_j = vec![0.011, 0.009, 0.012, 0.010, 0.008, 0.013, 0.007, 0.0095];
        (pos_i, pos_j, w_ij, vol_j)
    }

    /// Disordered 3-D stencil with 9 neighbours.
    fn disordered_stencil_3d() -> (Vec<f64>, Vec<Vec<f64>>, Vec<f64>, Vec<f64>) {
        let pos_i = vec![0.1, 0.2, -0.3];
        let offsets = [
            [0.10, 0.02, -0.03],
            [-0.08, 0.05, 0.04],
            [0.03, -0.11, 0.06],
            [-0.04, 0.09, -0.07],
            [0.07, 0.08, 0.09],
            [-0.11, -0.06, -0.05],
            [0.05, -0.04, -0.10],
            [-0.02, -0.03, 0.12],
            [0.09, -0.09, 0.01],
        ];
        let pos_j = offsets
            .iter()
            .map(|o| vec![pos_i[0] + o[0], pos_i[1] + o[1], pos_i[2] + o[2]])
            .collect();
        let w_ij = vec![0.9, 0.7, 0.6, 0.8, 0.5, 0.4, 1.1, 1.3, 0.75];
        let vol_j = vec![
            0.011, 0.009, 0.012, 0.010, 0.008, 0.013, 0.007, 0.0095, 0.0105,
        ];
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
    fn test_mls_correction_collinear_neighbours_invalid() {
        // All neighbours on the x-axis: the y direction is not spanned.
        let pos_i = vec![0.0, 0.0];
        let pos_j = vec![vec![0.1, 0.0], vec![-0.1, 0.0], vec![0.2, 0.0]];
        let w_ij = vec![1.0, 1.0, 0.5];
        let mls = mls_correction_matrix(&pos_i, &pos_j, &w_ij, 2);
        assert!(!mls.is_valid);
        assert!((mls.normalization - 2.5).abs() < 1e-14);
    }

    #[test]
    fn test_mls_correction_valid_for_regular_2d_stencil() {
        let pos_i = vec![0.0, 0.0];
        let mut pos_j = Vec::new();
        for ix in -1_i32..=1 {
            for iy in -1_i32..=1 {
                if ix != 0 || iy != 0 {
                    pos_j.push(vec![f64::from(ix) * 0.1, f64::from(iy) * 0.1]);
                }
            }
        }
        let w_ij = vec![1.0; pos_j.len()];
        let mls = mls_correction_matrix(&pos_i, &pos_j, &w_ij, 2);
        assert!(mls.is_valid);
        // A = diag-like for the symmetric stencil: A[1][1] = 6 * 0.01 = 0.06,
        // so (A^-1)[1][1] = 1 / 0.06 and the cross terms vanish.
        assert!(
            (mls.get(0, 1) - 1.0 / 0.06).abs() < 1e-9,
            "B00 = {}",
            mls.get(0, 1)
        );
        assert!((mls.get(1, 2) - 1.0 / 0.06).abs() < 1e-9);
        assert!(mls.get(0, 2).abs() < 1e-9);
        assert!(mls.get(1, 1).abs() < 1e-9);
    }

    #[test]
    fn test_mls_correction_normalization_positive() {
        let pos_i = vec![0.0, 0.0];
        let pos_j = vec![vec![0.5, 0.0], vec![-0.5, 0.0], vec![0.0, 0.5]];
        let w_ij = vec![0.5, 0.5, 0.5];
        let mls = mls_correction_matrix(&pos_i, &pos_j, &w_ij, 2);
        assert!(mls.is_valid);
        assert!((mls.normalization - 1.5).abs() < 1e-14);
    }

    #[test]
    fn test_mls_correction_dim_preserved() {
        let pos_i = vec![0.0, 0.0];
        let pos_j = vec![vec![1.0, 0.0]];
        let w_ij = vec![1.0];
        let mls = mls_correction_matrix(&pos_i, &pos_j, &w_ij, 2);
        assert_eq!(mls.dim, 2);
        assert!(!mls.is_valid, "a single neighbour cannot span 2-D");
    }

    #[test]
    fn test_mls_correction_rows_are_inverse_gradient_rows() {
        // B must be rows 1.. of A^-1: (B A)[d][k] = delta(d + 1, k).
        let (pos_i, pos_j, w_ij, vol_j) = disordered_stencil_2d();
        let weights: Vec<f64> = w_ij.iter().zip(vol_j.iter()).map(|(w, v)| w * v).collect();
        let mls = mls_correction_matrix(&pos_i, &pos_j, &weights, 2);
        assert!(mls.is_valid);

        let p_size = 3;
        let mut a = vec![0.0f64; p_size * p_size];
        for (pj, &wv) in pos_j.iter().zip(weights.iter()) {
            let p = [1.0, pj[0] - pos_i[0], pj[1] - pos_i[1]];
            for r in 0..p_size {
                for c in 0..p_size {
                    a[r * p_size + c] += wv * p[r] * p[c];
                }
            }
        }
        for d in 0..2 {
            for k in 0..p_size {
                let mut entry = 0.0;
                for c in 0..p_size {
                    entry += mls.get(d, c) * a[c * p_size + k];
                }
                let expected = if k == d + 1 { 1.0 } else { 0.0 };
                assert!(
                    (entry - expected).abs() < 1e-9,
                    "(B A)[{d}][{k}] = {entry}, expected {expected}"
                );
            }
        }
    }

    #[test]
    fn test_mls_correction_length_mismatch_invalid() {
        let pos_i = vec![0.0, 0.0];
        let pos_j = vec![vec![0.1, 0.0], vec![0.0, 0.1]];
        let mls = mls_correction_matrix(&pos_i, &pos_j, &[1.0], 2);
        assert!(!mls.is_valid);
        let mls = mls_correction_matrix(&[0.0], &pos_j, &[1.0, 1.0], 2);
        assert!(!mls.is_valid);
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
        // No neighbours: the MLS system is undetermined, so neither the zeroth
        // nor the first order condition can be established.
        let pos_i = vec![0.0, 0.0];
        let (zeroth, first) = consistency_check(&pos_i, &[], &[], &[], 2, 1e-8);
        assert!(!zeroth);
        assert!(!first);
    }

    #[test]
    fn test_consistency_check_singular_stencil_reports_false() {
        // Collinear neighbours cannot reproduce a gradient in the y direction.
        let pos_i = vec![0.0, 0.0];
        let pos_j = vec![vec![0.1, 0.0], vec![-0.1, 0.0]];
        let w_ij = vec![0.9, 0.9];
        let vol_j = vec![0.1, 0.1];
        assert_eq!(
            consistency_check(&pos_i, &pos_j, &w_ij, &vol_j, 2, 1e-6),
            (false, false)
        );
    }

    #[test]
    fn test_consistency_check_disordered_2d_stencil() {
        let (pos_i, pos_j, w_ij, vol_j) = disordered_stencil_2d();
        assert_eq!(
            consistency_check(&pos_i, &pos_j, &w_ij, &vol_j, 2, 1e-10),
            (true, true)
        );
    }

    #[test]
    fn test_consistency_check_disordered_3d_stencil() {
        let (pos_i, pos_j, w_ij, vol_j) = disordered_stencil_3d();
        assert_eq!(
            consistency_check(&pos_i, &pos_j, &w_ij, &vol_j, 3, 1e-10),
            (true, true)
        );
    }

    #[test]
    fn test_consistency_check_rejects_malformed_positions() {
        let pos_i = vec![0.0, 0.0];
        let pos_j = vec![vec![0.1], vec![0.0, 0.1], vec![-0.1, -0.1]];
        let w_ij = vec![1.0; 3];
        let vol_j = vec![0.1; 3];
        assert_eq!(
            consistency_check(&pos_i, &pos_j, &w_ij, &vol_j, 2, 1e-6),
            (false, false)
        );
        assert_eq!(
            consistency_check(&[0.0], &pos_j, &w_ij, &vol_j, 2, 1e-6),
            (false, false)
        );
    }

    // ── mls_gradient ──────────────────────────────────────────────────────

    #[test]
    fn test_mls_gradient_no_neighbours_is_none() {
        let pos_i = vec![0.0, 0.0];
        assert!(mls_gradient(1.0, &[], &pos_i, &[], &[], &[], 2).is_none());
    }

    #[test]
    fn test_mls_gradient_single_neighbour_is_none() {
        let pos_i = vec![0.0, 0.0];
        let pos_j = vec![vec![0.1, 0.0]];
        let grad = mls_gradient(1.0, &[2.0], &pos_i, &pos_j, &[0.9], &[0.1], 2);
        assert!(
            grad.is_none(),
            "one neighbour cannot determine a 2-D gradient"
        );
    }

    #[test]
    fn test_mls_gradient_length_mismatch_is_none() {
        let pos_i = vec![0.0, 0.0];
        let pos_j = vec![vec![0.1, 0.0], vec![0.0, 0.1], vec![-0.1, -0.1]];
        let grad = mls_gradient(0.0, &[0.0; 3], &pos_i, &pos_j, &[1.0; 3], &[0.1; 2], 2);
        assert!(grad.is_none());
    }

    #[test]
    fn test_mls_gradient_constant_field_zero_grad() {
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
        let grad = mls_gradient(5.0, &phi_j, &pos_i, &pos_j, &w_ij, &vol_j, 2)
            .expect("symmetric stencil spans 2-D");
        assert_eq!(grad.len(), 2);
        for g in &grad {
            assert!(g.abs() < 1e-10, "gradient component {g} not near zero");
        }
    }

    #[test]
    fn test_mls_gradient_x_field_on_six_neighbours() {
        // Regression: the previous formula summed A^-1 over the neighbours and
        // returned ~ (-2.4e-5, -1.2e-3) instead of (1, 0) for phi = x.
        let pos_i = vec![0.0, 0.0];
        let pos_j = vec![
            vec![0.10, 0.00],
            vec![-0.10, 0.02],
            vec![0.00, 0.10],
            vec![0.03, -0.10],
            vec![0.08, 0.07],
            vec![-0.06, -0.05],
        ];
        let w_ij = vec![0.9, 0.7, 0.8, 0.6, 0.5, 0.4];
        let vol_j = vec![0.01; 6];
        let phi_j: Vec<f64> = pos_j.iter().map(|p| p[0]).collect();
        let grad = mls_gradient(0.0, &phi_j, &pos_i, &pos_j, &w_ij, &vol_j, 2)
            .expect("six neighbours span 2-D");
        assert!((grad[0] - 1.0).abs() < 1e-10, "d/dx = {}", grad[0]);
        assert!(grad[1].abs() < 1e-10, "d/dy = {}", grad[1]);
    }

    #[test]
    fn test_mls_gradient_exact_for_affine_field_2d() {
        let (pos_i, pos_j, w_ij, vol_j) = disordered_stencil_2d();
        let field = |p: &[f64]| 3.0 + 2.0 * p[0] - 5.0 * p[1];
        let phi_j: Vec<f64> = pos_j.iter().map(|p| field(p)).collect();
        let grad = mls_gradient(field(&pos_i), &phi_j, &pos_i, &pos_j, &w_ij, &vol_j, 2)
            .expect("disordered stencil spans 2-D");
        assert!((grad[0] - 2.0).abs() < 1e-10, "d/dx = {}", grad[0]);
        assert!((grad[1] + 5.0).abs() < 1e-10, "d/dy = {}", grad[1]);
    }

    #[test]
    fn test_mls_gradient_exact_for_affine_field_3d() {
        let (pos_i, pos_j, w_ij, vol_j) = disordered_stencil_3d();
        let field = |p: &[f64]| -1.0 + p[0] + 2.0 * p[1] - 3.0 * p[2];
        let phi_j: Vec<f64> = pos_j.iter().map(|p| field(p)).collect();
        let grad = mls_gradient(field(&pos_i), &phi_j, &pos_i, &pos_j, &w_ij, &vol_j, 3)
            .expect("disordered stencil spans 3-D");
        assert!((grad[0] - 1.0).abs() < 1e-10, "d/dx = {}", grad[0]);
        assert!((grad[1] - 2.0).abs() < 1e-10, "d/dy = {}", grad[1]);
        assert!((grad[2] + 3.0).abs() < 1e-10, "d/dz = {}", grad[2]);
    }

    #[test]
    fn test_mls_gradient_regular_stencil_matches_central_difference() {
        // For a symmetric 3x3 stencil with unit weights the MLS gradient of a
        // quadratic field phi = x^2 + y at the origin equals the exact (0, 1)
        // because the odd moments cancel.
        let pos_i = vec![0.0, 0.0];
        let mut pos_j = Vec::new();
        for ix in -1_i32..=1 {
            for iy in -1_i32..=1 {
                if ix != 0 || iy != 0 {
                    pos_j.push(vec![f64::from(ix) * 0.1, f64::from(iy) * 0.1]);
                }
            }
        }
        let w_ij = vec![1.0; pos_j.len()];
        let vol_j = vec![0.01; pos_j.len()];
        let phi_j: Vec<f64> = pos_j.iter().map(|p| p[0] * p[0] + p[1]).collect();
        let grad = mls_gradient(0.0, &phi_j, &pos_i, &pos_j, &w_ij, &vol_j, 2)
            .expect("3x3 stencil spans 2-D");
        assert!(grad[0].abs() < 1e-10, "d/dx = {}", grad[0]);
        assert!((grad[1] - 1.0).abs() < 1e-10, "d/dy = {}", grad[1]);
    }

    #[test]
    fn test_mls_gradient_is_length_unit_independent() {
        // The same stencil scaled to nano-metre size (positions x 1e-9, kernel
        // values x 1e27, volumes x 1e-27) must still be recognised as regular
        // and reproduce phi = x with unit slope.
        let scale = 1e-9;
        let (pos_i, pos_j, w_ij, vol_j) = disordered_stencil_2d();
        let pos_i_s: Vec<f64> = pos_i.iter().map(|x| x * scale).collect();
        let pos_j_s: Vec<Vec<f64>> = pos_j
            .iter()
            .map(|p| p.iter().map(|x| x * scale).collect())
            .collect();
        let w_s: Vec<f64> = w_ij.iter().map(|w| w * 1e27).collect();
        let v_s: Vec<f64> = vol_j.iter().map(|v| v * 1e-27).collect();
        let phi_j: Vec<f64> = pos_j_s.iter().map(|p| p[0]).collect();
        let grad = mls_gradient(pos_i_s[0], &phi_j, &pos_i_s, &pos_j_s, &w_s, &v_s, 2)
            .expect("nano-scale stencil is still regular");
        assert!((grad[0] - 1.0).abs() < 1e-9, "d/dx = {}", grad[0]);
        assert!(grad[1].abs() < 1e-9, "d/dy = {}", grad[1]);
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
