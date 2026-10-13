// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Spectral analysis of operators represented as finite matrices.
//!
//! Power iteration for the dominant eigenvalue, shift-and-invert iteration for
//! the eigenvalue nearest a given shift, and the spectral-norm resolvent.

// ─────────────────────────────────────────────────────────────────────────────
// OperatorSpectrum
// ─────────────────────────────────────────────────────────────────────────────

/// Spectral analysis of compact operators represented as finite matrices.
///
/// Provides eigenvalue computation via power iteration, spectral radius
/// estimation, and resolvent norm.
#[derive(Debug, Clone)]
pub struct OperatorSpectrum {
    /// Matrix representation of the operator (row-major, square).
    pub matrix: Vec<Vec<f64>>,
    /// Size of the square matrix.
    pub n: usize,
}

impl OperatorSpectrum {
    /// Create an `OperatorSpectrum` from a square matrix.
    ///
    /// Panics if the matrix is not square.
    pub fn new(matrix: Vec<Vec<f64>>) -> Self {
        let n = matrix.len();
        assert!(
            matrix.iter().all(|row| row.len() == n),
            "matrix must be square"
        );
        Self { matrix, n }
    }

    /// Estimate the largest eigenvalue (spectral radius) via the power method.
    ///
    /// Returns the dominant eigenvalue after `max_iter` iterations.
    pub fn spectral_radius(&self, max_iter: usize) -> f64 {
        if self.n == 0 {
            return 0.0;
        }
        let mut v = vec![1.0 / (self.n as f64).sqrt(); self.n];
        let mut lambda = 0.0_f64;
        for _ in 0..max_iter {
            let w = self.matvec(&v);
            let norm: f64 = w.iter().map(|x| x * x).sum::<f64>().sqrt();
            if norm < 1e-15 {
                break;
            }
            lambda = w.iter().zip(v.iter()).map(|(wi, vi)| wi * vi).sum::<f64>();
            for (vi, wi) in v.iter_mut().zip(w.iter()) {
                *vi = wi / norm;
            }
        }
        lambda
    }

    /// Estimate the eigenvalue of the operator closest to `shift` by
    /// shift-and-invert (inverse) power iteration.
    ///
    /// Each iteration solves `(A − shift·I) w = v` with a single LU
    /// factorisation of the shifted matrix and normalises `w`.  The iterates
    /// converge to the eigenvector whose eigenvalue is nearest to `shift`
    /// (the eigenvalue with the largest modulus of `(A − shift·I)⁻¹`), and the
    /// returned value is the Rayleigh quotient `vᵀ A v` of the final iterate.
    /// Iteration stops after `max_iter` steps or once the eigen-residual
    /// `‖A v − λ v‖` drops below `1e-12` relative to the matrix scale.
    ///
    /// * With `shift = 0.0` this is the eigenvalue of smallest *magnitude*; it is
    ///   the algebraically smallest eigenvalue only when the spectrum is
    ///   non-negative.  Any `shift` below the spectrum yields the smallest
    ///   eigenvalue, any `shift` above it the largest.
    /// * The start vector is a fixed, non-symmetric pseudo-random vector so that
    ///   operators with structured eigenvectors (for example circulant matrices,
    ///   whose all-ones mode would otherwise trap the iteration) are handled.
    /// * If `shift` coincides with an eigenvalue the shifted matrix is singular
    ///   and `shift` itself is returned.
    /// * Real matrices whose nearest eigenvalues form a complex-conjugate pair
    ///   (equal distance from `shift`) do not converge; the estimate after
    ///   `max_iter` steps is returned.
    ///
    /// Returns `0.0` for an empty operator and `NaN` if the matrix or the shift
    /// contains a non-finite value.
    pub fn smallest_eigenvalue(&self, shift: f64, max_iter: usize) -> f64 {
        shift_invert_eigenvalue(&self.matrix, shift, max_iter)
    }

    /// Compute the trace of the matrix (sum of diagonal elements).
    pub fn trace(&self) -> f64 {
        (0..self.n).map(|i| self.matrix[i][i]).sum()
    }

    /// Compute the Frobenius norm `‖A‖_F = √(Σ aᵢⱼ²)`.
    pub fn frobenius_norm(&self) -> f64 {
        self.matrix
            .iter()
            .flat_map(|row| row.iter())
            .map(|x| x * x)
            .sum::<f64>()
            .sqrt()
    }

    /// Compute the spectral-norm resolvent `‖(A − λI)⁻¹‖₂` at the point `lambda`.
    ///
    /// For the Euclidean operator norm this equals `1 / σ_min(A − λI)`, where
    /// `σ_min` is the smallest singular value of the shifted matrix, evaluated
    /// from a full singular value decomposition.  For a normal matrix it is
    /// `1 / dist(λ, spectrum)`; for a non-normal matrix it can be much larger.
    ///
    /// Returns [`f64::INFINITY`] when `λ` is an eigenvalue (that is, when
    /// `σ_min` is not larger than `1e-14 · σ_max`), `0.0` for an empty operator
    /// and `NaN` if the matrix or `lambda` contains a non-finite value.
    pub fn resolvent_norm(&self, lambda: f64) -> f64 {
        if self.n == 0 {
            return 0.0;
        }
        if !lambda.is_finite() || self.matrix.iter().flatten().any(|a| !a.is_finite()) {
            return f64::NAN;
        }
        let n = self.n;
        let shifted = nalgebra::DMatrix::<f64>::from_fn(n, n, |i, j| {
            if i == j {
                self.matrix[i][j] - lambda
            } else {
                self.matrix[i][j]
            }
        });
        let Some(svd) = nalgebra::SVD::try_new(shifted, false, false, f64::EPSILON, 1000 * n)
        else {
            return f64::NAN;
        };
        let sigma_max = svd
            .singular_values
            .iter()
            .fold(0.0_f64, |acc, &sv| acc.max(sv));
        let sigma_min = svd
            .singular_values
            .iter()
            .fold(f64::INFINITY, |acc, &sv| acc.min(sv));
        if sigma_min <= 1e-14 * sigma_max {
            return f64::INFINITY;
        }
        1.0 / sigma_min
    }

    /// Apply the operator (matrix-vector multiplication) to `v`.
    pub fn apply(&self, v: &[f64]) -> Vec<f64> {
        self.matvec(v)
    }

    /// Compute all diagonal elements (approximate eigenvalues of a diagonal operator).
    pub fn diagonal(&self) -> Vec<f64> {
        (0..self.n).map(|i| self.matrix[i][i]).collect()
    }

    fn matvec(&self, v: &[f64]) -> Vec<f64> {
        self.matrix
            .iter()
            .map(|row| row.iter().zip(v.iter()).map(|(a, x)| a * x).sum())
            .collect()
    }
}

/// Shift-and-invert iteration for the eigenvalue of `matrix` closest to `shift`.
///
/// See [`OperatorSpectrum::smallest_eigenvalue`] for the contract.  The matrix is
/// scaled by its largest entry so that the absolute pivot tolerance of
/// [`crate::linalg::lu_factor_n`] is a relative one, and the shifted matrix is
/// factorised once for all iterations.
fn shift_invert_eigenvalue(matrix: &[Vec<f64>], shift: f64, max_iter: usize) -> f64 {
    let n = matrix.len();
    if n == 0 {
        return 0.0;
    }
    if !shift.is_finite() || matrix.iter().flatten().any(|a| !a.is_finite()) {
        return f64::NAN;
    }
    let scale = matrix
        .iter()
        .flatten()
        .fold(0.0_f64, |acc, a| acc.max(a.abs()));
    if scale == 0.0 {
        // The zero operator: every eigenvalue is 0.
        return 0.0;
    }
    let scaled: Vec<Vec<f64>> = matrix
        .iter()
        .map(|row| row.iter().map(|a| a / scale).collect())
        .collect();
    let sigma = shift / scale;
    let mut shifted = scaled.clone();
    for (i, row) in shifted.iter_mut().enumerate() {
        row[i] -= sigma;
    }
    let Some((lu, perm)) = crate::linalg::lu_factor_n(&shifted) else {
        // Singular shifted matrix: `shift` is an eigenvalue.
        return shift;
    };
    // Fixed low-discrepancy start vector (golden-ratio sequence).  It is not
    // symmetric and has a non-zero component along every mode of structured
    // operators such as circulant matrices.
    let mut v: Vec<f64> = (0..n)
        .map(|i| 0.1 + ((i + 1) as f64 * 0.618_033_988_749_894_9).fract())
        .collect();
    normalise_in_place(&mut v);
    let mut lambda = sigma;
    for _ in 0..max_iter {
        let permuted: Vec<f64> = perm.iter().map(|&p| v[p]).collect();
        let forward = crate::linalg::forward_sub_n(&lu, &permuted);
        let Some(mut w) = crate::linalg::backward_sub_n(&lu, &forward) else {
            return shift;
        };
        if !normalise_in_place(&mut w) {
            break;
        }
        v = w;
        let av = crate::linalg::matvec_n(&scaled, &v);
        lambda = v.iter().zip(av.iter()).map(|(vi, avi)| vi * avi).sum();
        let residual = av
            .iter()
            .zip(v.iter())
            .map(|(avi, vi)| (avi - lambda * vi).powi(2))
            .sum::<f64>()
            .sqrt();
        if residual <= 1e-12 {
            break;
        }
    }
    lambda * scale
}

/// Scale `v` to unit Euclidean length; returns `false` (leaving `v` untouched)
/// if its norm is zero or not finite.
fn normalise_in_place(v: &mut [f64]) -> bool {
    let norm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
    if !norm.is_finite() || norm == 0.0 {
        return false;
    }
    for x in v.iter_mut() {
        *x /= norm;
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── OperatorSpectrum ──────────────────────────────────────────────────────

    #[test]
    fn test_operator_spectrum_trace() {
        let m = vec![vec![2.0, 1.0], vec![0.0, 3.0]];
        let op = OperatorSpectrum::new(m);
        assert!((op.trace() - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_operator_spectrum_frobenius() {
        let m = vec![vec![3.0, 4.0], vec![0.0, 0.0]];
        let op = OperatorSpectrum::new(m);
        assert!((op.frobenius_norm() - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_operator_spectrum_apply() {
        let m = vec![vec![2.0, 0.0], vec![0.0, 3.0]];
        let op = OperatorSpectrum::new(m);
        let v = vec![1.0, 2.0];
        let w = op.apply(&v);
        assert!((w[0] - 2.0).abs() < 1e-12);
        assert!((w[1] - 6.0).abs() < 1e-12);
    }

    #[test]
    fn test_operator_spectrum_spectral_radius() {
        // Diagonal matrix: spectral radius = max eigenvalue
        let m = vec![vec![5.0, 0.0], vec![0.0, 2.0]];
        let op = OperatorSpectrum::new(m);
        let rho = op.spectral_radius(200);
        assert!((rho - 5.0).abs() < 0.1, "spectral radius = {rho}");
    }

    #[test]
    fn test_operator_spectrum_diagonal() {
        let m = vec![vec![7.0, 0.0], vec![0.0, 3.0]];
        let op = OperatorSpectrum::new(m);
        let d = op.diagonal();
        assert!((d[0] - 7.0).abs() < 1e-12);
        assert!((d[1] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_operator_spectrum_resolvent_finite() {
        let m = vec![vec![1.0, 0.0], vec![0.0, 2.0]];
        let op = OperatorSpectrum::new(m);
        // At a regular point (e.g., λ=10) resolvent should be finite
        let r = op.resolvent_norm(10.0);
        assert!(r.is_finite() && r > 0.0, "resolvent = {r}");
    }

    #[test]
    fn test_smallest_eigenvalue_returns_eigenvalue_nearest_zero() {
        // The forward power method would return 10 here.
        let op = OperatorSpectrum::new(vec![vec![1.0, 0.0], vec![0.0, 10.0]]);
        let lam = op.smallest_eigenvalue(0.0, 200);
        assert!((lam - 1.0).abs() < 1e-10, "smallest eigenvalue = {lam}");
        // The dominant eigenvalue is still available from the power method.
        assert!((op.spectral_radius(200) - 10.0).abs() < 1e-6);
    }

    #[test]
    fn test_smallest_eigenvalue_shift_selects_nearest_eigenvalue() {
        let op = OperatorSpectrum::new(vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 4.0, 0.0],
            vec![0.0, 0.0, 10.0],
        ]);
        assert!((op.smallest_eigenvalue(0.5, 200) - 1.0).abs() < 1e-10);
        assert!((op.smallest_eigenvalue(3.0, 200) - 4.0).abs() < 1e-10);
        assert!((op.smallest_eigenvalue(9.0, 200) - 10.0).abs() < 1e-10);
        // A shift that is itself an eigenvalue is returned unchanged.
        assert_eq!(op.smallest_eigenvalue(4.0, 200), 4.0);
    }

    #[test]
    fn test_smallest_eigenvalue_symmetric_matrix_with_ones_eigenvector() {
        // Eigenvalues 1 (v = [1,-1]) and 3 (v = [1,1]).  A start vector equal to
        // the all-ones vector would converge to 3.
        let op = OperatorSpectrum::new(vec![vec![2.0, 1.0], vec![1.0, 2.0]]);
        let lam = op.smallest_eigenvalue(0.0, 200);
        assert!((lam - 1.0).abs() < 1e-10, "got {lam}");
    }

    #[test]
    fn test_smallest_eigenvalue_circulant_matrix() {
        // Periodic 1-D Laplacian on 6 nodes: eigenvalues 2 - 2cos(2πk/6) =
        // {0, 1, 1, 3, 3, 4}.  The all-ones vector is the eigenvector of 0.
        let n = 6;
        let mut m = vec![vec![0.0; n]; n];
        for i in 0..n {
            m[i][i] = 2.0;
            m[i][(i + 1) % n] = -1.0;
            m[i][(i + n - 1) % n] = -1.0;
        }
        let op = OperatorSpectrum::new(m);
        let lam = op.smallest_eigenvalue(1.2, 500);
        assert!(
            (lam - 1.0).abs() < 1e-9,
            "closest to 1.2 should be 1, got {lam}"
        );
        let lam = op.smallest_eigenvalue(3.4, 500);
        assert!(
            (lam - 3.0).abs() < 1e-9,
            "closest to 3.4 should be 3, got {lam}"
        );
        let lam = op.smallest_eigenvalue(3.9, 500);
        assert!(
            (lam - 4.0).abs() < 1e-9,
            "closest to 3.9 should be 4, got {lam}"
        );
    }

    #[test]
    fn test_smallest_eigenvalue_non_symmetric_matrix() {
        // Upper triangular: eigenvalues 2 and 3 with non-orthogonal eigenvectors.
        let op = OperatorSpectrum::new(vec![vec![2.0, 5.0], vec![0.0, 3.0]]);
        assert!((op.smallest_eigenvalue(0.0, 200) - 2.0).abs() < 1e-9);
        assert!((op.smallest_eigenvalue(2.9, 200) - 3.0).abs() < 1e-9);
    }

    #[test]
    fn test_smallest_eigenvalue_is_scale_invariant() {
        for scale in [1e-20, 1.0, 1e20] {
            let op = OperatorSpectrum::new(vec![vec![scale, 0.0], vec![0.0, 10.0 * scale]]);
            let lam = op.smallest_eigenvalue(0.0, 200);
            assert!(
                (lam / scale - 1.0).abs() < 1e-9,
                "scale {scale}: eigenvalue {lam}"
            );
        }
    }

    #[test]
    fn test_smallest_eigenvalue_degenerate_inputs() {
        assert_eq!(
            OperatorSpectrum::new(vec![]).smallest_eigenvalue(0.0, 10),
            0.0
        );
        let one = OperatorSpectrum::new(vec![vec![-3.5]]);
        assert!((one.smallest_eigenvalue(0.0, 10) + 3.5).abs() < 1e-12);
        let zero = OperatorSpectrum::new(vec![vec![0.0, 0.0], vec![0.0, 0.0]]);
        assert_eq!(zero.smallest_eigenvalue(0.0, 10), 0.0);
        let bad = OperatorSpectrum::new(vec![vec![f64::NAN]]);
        assert!(bad.smallest_eigenvalue(0.0, 10).is_nan());
        let op = OperatorSpectrum::new(vec![vec![1.0, 0.0], vec![0.0, 10.0]]);
        assert!(op.smallest_eigenvalue(f64::INFINITY, 10).is_nan());
    }

    #[test]
    fn test_resolvent_norm_equals_inverse_distance_to_spectrum() {
        // Normal (diagonal) operator: ||(A - λ)⁻¹|| = 1 / dist(λ, {1, 2}).
        let op = OperatorSpectrum::new(vec![vec![1.0, 0.0], vec![0.0, 2.0]]);
        let r = op.resolvent_norm(1.001);
        assert!((r - 1000.0).abs() < 1e-6, "resolvent = {r}");
        let r = op.resolvent_norm(10.0);
        assert!((r - 0.125).abs() < 1e-12, "resolvent = {r}");
        // The resolvent blows up as λ approaches an eigenvalue.
        assert!(op.resolvent_norm(1.0 + 1e-8) > op.resolvent_norm(1.0 + 1e-4));
    }

    #[test]
    fn test_resolvent_norm_infinite_at_eigenvalue_and_zero_when_empty() {
        let op = OperatorSpectrum::new(vec![vec![1.0, 0.0], vec![0.0, 2.0]]);
        assert!(op.resolvent_norm(1.0).is_infinite());
        assert!(op.resolvent_norm(2.0).is_infinite());
        assert_eq!(OperatorSpectrum::new(vec![]).resolvent_norm(1.0), 0.0);
        let bad = OperatorSpectrum::new(vec![vec![f64::INFINITY]]);
        assert!(bad.resolvent_norm(0.0).is_nan());
        assert!(op.resolvent_norm(f64::NAN).is_nan());
    }

    #[test]
    fn test_resolvent_norm_non_normal_operator_exceeds_spectral_estimate() {
        // A = [[1, 100], [0, 1]] has the single eigenvalue 1, so at λ = 0 the
        // spectral estimate 1/dist(λ, σ(A)) is 1, but the true resolvent norm
        // ||A⁻¹||₂ = 1/σ_min(A) is about 100.
        let op = OperatorSpectrum::new(vec![vec![1.0, 100.0], vec![0.0, 1.0]]);
        let r = op.resolvent_norm(0.0);
        // A⁻¹ = [[1, -100], [0, 1]]; its largest singular value solves
        // σ² = (2 + 100² + sqrt((2 + 100²)² - 4)) / 2.
        let t = 2.0 + 100.0_f64.powi(2);
        let expected = ((t + (t * t - 4.0).sqrt()) / 2.0).sqrt();
        assert!((r - expected).abs() < 1e-8 * expected, "{r} vs {expected}");
        assert!(r > 99.0);
    }
}
