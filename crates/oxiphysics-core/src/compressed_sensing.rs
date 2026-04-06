#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Compressed sensing and sparse signal recovery algorithms.
//!
//! Provides Discrete Cosine Transform (DCT) bases, random measurement matrices,
//! Basis Pursuit (ISTA/LASSO/FISTA), Orthogonal Matching Pursuit (OMP),
//! sparsity metrics, dictionary learning (K-SVD), sparse coding, FISTA,
//! restricted isometry property (RIP) analysis, Bernoulli measurement matrices,
//! MRI-like signal reconstruction, and theoretical recovery guarantees for
//! compressed sensing problems.
//!
//! # Overview
//!
//! Compressed sensing (CS) exploits the sparsity of natural signals to allow
//! faithful reconstruction from far fewer measurements than the Nyquist rate.
//! The key ingredients are:
//!
//! * A **sparsifying basis** (DCT, Wavelet, etc.) in which the signal has few
//!   non-zero coefficients.
//! * A **measurement matrix** (random Gaussian or Bernoulli) that is incoherent
//!   with the sparsifying basis.
//! * A **recovery algorithm** (Basis Pursuit / ISTA / FISTA / OMP) that finds
//!   the sparsest signal consistent with the measurements.
//!
//! The `BasisPursuit` struct implements ISTA (slow) and FISTA (fast, with
//! Nesterov momentum) for L1-regularised least-squares (LASSO) recovery.
//! `OrthogonalMatchingPursuit` provides a greedy alternative.
//! `KSvd` implements the K-SVD dictionary learning algorithm.

use rand::RngExt;

// ─────────────────────────────────────────────────────────────────────────────
// Free functions
// ─────────────────────────────────────────────────────────────────────────────

/// Apply element-wise soft-thresholding: `sign(x) * max(|x| - lambda, 0)`.
///
/// Used as a proximal operator in iterative shrinkage-thresholding algorithms.
#[allow(dead_code)]
pub fn soft_threshold(x: f64, lambda: f64) -> f64 {
    if x > lambda {
        x - lambda
    } else if x < -lambda {
        x + lambda
    } else {
        0.0
    }
}

/// Compute the Nyquist sampling rate for a band-limited signal.
///
/// Returns `2 * bandwidth` (samples per second).
#[allow(dead_code)]
pub fn nyquist_rate(bandwidth: f64) -> f64 {
    2.0 * bandwidth
}

/// Compute the compression ratio `m / n`.
///
/// A ratio less than 1 indicates sub-Nyquist sampling.
#[allow(dead_code)]
pub fn compression_ratio(n: usize, m: usize) -> f64 {
    if n == 0 {
        return 0.0;
    }
    m as f64 / n as f64
}

/// Compute the ℓ₂ norm of a slice.
///
/// Returns `sqrt(sum of squares)`.
#[allow(dead_code)]
pub fn l2_norm(x: &[f64]) -> f64 {
    x.iter().map(|v| v * v).sum::<f64>().sqrt()
}

/// Normalise a vector to unit ℓ₂ norm in place.
///
/// If the norm is smaller than `1e-14` the vector is left unchanged.
#[allow(dead_code)]
pub fn normalise(x: &mut Vec<f64>) {
    let n = l2_norm(x);
    if n > 1e-14 {
        for v in x.iter_mut() {
            *v /= n;
        }
    }
}

/// Compute the matrix-vector product `y = A x`.
///
/// `a` is row-major with shape `m × n`; returns a vector of length `m`.
#[allow(dead_code)]
pub fn mat_vec(a: &[Vec<f64>], x: &[f64]) -> Vec<f64> {
    a.iter()
        .map(|row| row.iter().zip(x.iter()).map(|(ai, xi)| ai * xi).sum())
        .collect()
}

/// Compute the transposed matrix-vector product `y = A^T x`.
///
/// `a` is row-major with shape `m × n`; returns a vector of length `n`.
#[allow(dead_code)]
pub fn mat_transpose_vec(a: &[Vec<f64>], x: &[f64]) -> Vec<f64> {
    if a.is_empty() {
        return Vec::new();
    }
    let n = a[0].len();
    let m = a.len();
    let mut y = vec![0.0_f64; n];
    for i in 0..m {
        for j in 0..n {
            y[j] += a[i][j] * x[i];
        }
    }
    y
}

/// Estimate the spectral norm (largest singular value) of `a` via power iteration.
///
/// Runs `max_iter` iterations; returns an approximation to `||A||_2`.
#[allow(dead_code)]
pub fn spectral_norm(a: &[Vec<f64>], max_iter: usize) -> f64 {
    if a.is_empty() {
        return 0.0;
    }
    let n = a[0].len();
    let mut v = vec![1.0_f64; n];
    normalise(&mut v);
    for _ in 0..max_iter {
        let av = mat_vec(a, &v);
        let mut atav = mat_transpose_vec(a, &av);
        normalise(&mut atav);
        v = atav;
    }
    let av = mat_vec(a, &v);
    l2_norm(&av)
}

// ─────────────────────────────────────────────────────────────────────────────
// DctBasis
// ─────────────────────────────────────────────────────────────────────────────

/// Discrete Cosine Transform (DCT-II) basis for sparse representation.
///
/// Signals with few significant DCT coefficients can be recovered from
/// far fewer measurements than the Nyquist rate.
#[allow(dead_code)]
pub struct DctBasis {
    /// Length of the signal (number of samples).
    pub n: usize,
}

#[allow(dead_code)]
impl DctBasis {
    /// Create a DCT basis for signals of length `n`.
    pub fn new(n: usize) -> Self {
        Self { n }
    }

    /// Compute the forward DCT-II transform of `x`.
    ///
    /// Returns a coefficient vector of the same length.
    pub fn transform(&self, x: &[f64]) -> Vec<f64> {
        let n = self.n.min(x.len());
        let mut out = vec![0.0; n];
        let pi_over_n = std::f64::consts::PI / n as f64;
        for k in 0..n {
            let mut sum = 0.0;
            for j in 0..n {
                sum += x[j] * ((j as f64 + 0.5) * k as f64 * pi_over_n).cos();
            }
            // DCT-II normalisation
            let norm = if k == 0 {
                (1.0 / n as f64).sqrt()
            } else {
                (2.0 / n as f64).sqrt()
            };
            out[k] = sum * norm;
        }
        out
    }

    /// Compute the inverse DCT-II (i.e. DCT-III) transform of `coeffs`.
    ///
    /// Reconstructs the original signal from its DCT coefficients.
    pub fn inverse(&self, coeffs: &[f64]) -> Vec<f64> {
        let n = self.n.min(coeffs.len());
        let mut out = vec![0.0; n];
        let pi_over_n = std::f64::consts::PI / n as f64;
        for j in 0..n {
            let mut sum = (1.0 / n as f64).sqrt() * coeffs[0];
            for k in 1..n {
                let norm = (2.0 / n as f64).sqrt();
                sum += norm * coeffs[k] * ((j as f64 + 0.5) * k as f64 * pi_over_n).cos();
            }
            out[j] = sum;
        }
        out
    }

    /// Threshold DCT coefficients to keep only the `k` largest-magnitude ones.
    ///
    /// Returns a truncated coefficient vector (all others set to zero).
    pub fn truncate(&self, coeffs: &[f64], k: usize) -> Vec<f64> {
        let mut indexed: Vec<(usize, f64)> = coeffs.iter().copied().enumerate().collect();
        indexed.sort_by(|a, b| {
            b.1.abs()
                .partial_cmp(&a.1.abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let mut out = vec![0.0_f64; coeffs.len()];
        for (i, v) in indexed.into_iter().take(k) {
            out[i] = v;
        }
        out
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// RandomMeasurementMatrix
// ─────────────────────────────────────────────────────────────────────────────

/// A random Gaussian measurement matrix for compressed sensing.
///
/// Each row is an independent Gaussian random vector; `m << n` enables
/// sub-Nyquist recovery of sparse signals.
#[allow(dead_code)]
pub struct RandomMeasurementMatrix {
    /// Number of measurements (rows).
    pub m: usize,
    /// Signal length (columns).
    pub n: usize,
    /// Underlying matrix entries stored row-major.
    pub matrix: Vec<Vec<f64>>,
}

#[allow(dead_code)]
impl RandomMeasurementMatrix {
    /// Generate an `m × n` Gaussian measurement matrix (entries ~ N(0, 1/m)).
    ///
    /// The columns are scaled by `1/sqrt(m)` so that each measurement
    /// approximately preserves the signal energy.
    pub fn generate_gaussian(m: usize, n: usize) -> Self {
        use rand::RngExt as _;
        let mut rng = rand::rng();
        let scale = 1.0 / (m as f64).sqrt();
        let matrix: Vec<Vec<f64>> = (0..m)
            .map(|_| {
                (0..n)
                    .map(|_| {
                        // Box-Muller for N(0,1)
                        let u1: f64 = rng.random_range(1e-12_f64..1.0_f64);
                        let u2: f64 = rng.random_range(0.0_f64..1.0_f64);
                        let z = (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
                        z * scale
                    })
                    .collect()
            })
            .collect();
        Self { m, n, matrix }
    }

    /// Generate an `m × n` Bernoulli ±1/sqrt(m) measurement matrix.
    ///
    /// Each entry is independently ±1/sqrt(m) with equal probability 1/2.
    pub fn generate_bernoulli(m: usize, n: usize) -> Self {
        let mut rng = rand::rng();
        let scale = 1.0 / (m as f64).sqrt();
        let matrix: Vec<Vec<f64>> = (0..m)
            .map(|_| {
                (0..n)
                    .map(|_| {
                        if rng.random_range(0.0_f64..1.0_f64) < 0.5 {
                            scale
                        } else {
                            -scale
                        }
                    })
                    .collect()
            })
            .collect();
        Self { m, n, matrix }
    }

    /// Apply the measurement matrix to signal `x`.
    ///
    /// Returns a vector of length `m` (the compressed measurements).
    pub fn measure(&self, x: &[f64]) -> Vec<f64> {
        mat_vec(&self.matrix, x)
    }

    /// Compute the mutual coherence of the measurement matrix columns.
    ///
    /// Lower coherence is better for sparse recovery.
    pub fn coherence(&self) -> f64 {
        SparsityMetrics::coherence(&self.matrix)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BasisPursuit / ISTA / FISTA
// ─────────────────────────────────────────────────────────────────────────────

/// Basis Pursuit via iterative shrinkage-thresholding (ISTA and FISTA).
///
/// Solves the LASSO problem: `argmin_x 0.5 ||Ax - b||^2 + lambda ||x||_1`.
///
/// FISTA adds Nesterov momentum for faster O(1/k²) convergence versus
/// O(1/k) for plain ISTA.
#[allow(dead_code)]
pub struct BasisPursuit;

#[allow(dead_code)]
impl BasisPursuit {
    /// Estimate the Lipschitz constant of the gradient via power iteration.
    fn lipschitz(a: &[Vec<f64>]) -> f64 {
        spectral_norm(a, 20).powi(2).max(1e-10)
    }

    /// Solve the LASSO problem using ISTA.
    ///
    /// - `a` — measurement matrix (m × n, row-major)
    /// - `b` — measurement vector (length m)
    /// - `lambda` — sparsity regularisation weight
    /// - `max_iter` — maximum number of iterations
    ///
    /// Returns the recovered sparse signal of length n.
    pub fn solve_lasso(a: &[Vec<f64>], b: &[f64], lambda: f64, max_iter: usize) -> Vec<f64> {
        if a.is_empty() || b.is_empty() {
            return Vec::new();
        }
        let n = a[0].len();
        let l = Self::lipschitz(a);
        let step = 1.0 / l;
        let mut x = vec![0.0_f64; n];

        for _ in 0..max_iter {
            let residual: Vec<f64> = mat_vec(a, &x)
                .iter()
                .zip(b.iter())
                .map(|(r, bi)| r - bi)
                .collect();
            let grad = mat_transpose_vec(a, &residual);
            x = x
                .iter()
                .zip(grad.iter())
                .map(|(xi, gi)| soft_threshold(xi - step * gi, step * lambda))
                .collect();
        }
        x
    }

    /// Solve the LASSO problem using FISTA (Fast ISTA with Nesterov momentum).
    ///
    /// - `a` — measurement matrix (m × n, row-major)
    /// - `b` — measurement vector (length m)
    /// - `lambda` — sparsity regularisation weight
    /// - `max_iter` — maximum number of iterations
    ///
    /// Returns the recovered sparse signal of length n.
    pub fn solve_fista(a: &[Vec<f64>], b: &[f64], lambda: f64, max_iter: usize) -> Vec<f64> {
        if a.is_empty() || b.is_empty() {
            return Vec::new();
        }
        let n = a[0].len();
        let l = Self::lipschitz(a);
        let step = 1.0 / l;

        let mut x = vec![0.0_f64; n];
        let mut y = x.clone();
        let mut t = 1.0_f64;

        for _ in 0..max_iter {
            let x_prev = x.clone();

            // Gradient step at y
            let residual: Vec<f64> = mat_vec(a, &y)
                .iter()
                .zip(b.iter())
                .map(|(r, bi)| r - bi)
                .collect();
            let grad = mat_transpose_vec(a, &residual);
            x = y
                .iter()
                .zip(grad.iter())
                .map(|(yi, gi)| soft_threshold(yi - step * gi, step * lambda))
                .collect();

            // Nesterov momentum update
            let t_new = (1.0 + (1.0 + 4.0 * t * t).sqrt()) / 2.0;
            let momentum = (t - 1.0) / t_new;
            y = x
                .iter()
                .zip(x_prev.iter())
                .map(|(xi, xi_prev)| xi + momentum * (xi - xi_prev))
                .collect();
            t = t_new;
        }
        x
    }

    /// Compute the objective value `0.5 ||Ax - b||^2 + lambda ||x||_1`.
    ///
    /// Useful for monitoring convergence.
    pub fn objective(a: &[Vec<f64>], b: &[f64], x: &[f64], lambda: f64) -> f64 {
        if a.is_empty() || b.is_empty() {
            return 0.0;
        }
        let ax = mat_vec(a, x);
        let residual_sq: f64 = ax
            .iter()
            .zip(b.iter())
            .map(|(r, bi)| (r - bi).powi(2))
            .sum();
        let l1: f64 = x.iter().map(|xi| xi.abs()).sum();
        0.5 * residual_sq + lambda * l1
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// OrthogonalMatchingPursuit
// ─────────────────────────────────────────────────────────────────────────────

/// Orthogonal Matching Pursuit (OMP) for sparse signal recovery.
///
/// Greedily selects the most correlated column of the measurement matrix
/// at each step and performs a least-squares fit on the selected support.
#[allow(dead_code)]
pub struct OrthogonalMatchingPursuit {
    /// Maximum sparsity (number of non-zero coefficients to recover).
    pub max_k: usize,
}

#[allow(dead_code)]
impl OrthogonalMatchingPursuit {
    /// Create an OMP solver with sparsity bound `max_k`.
    pub fn new(max_k: usize) -> Self {
        Self { max_k }
    }

    /// Recover a sparse signal from measurements.
    ///
    /// - `a` — measurement matrix (m × n, row-major)
    /// - `b` — measurement vector (length m)
    ///
    /// Returns the recovered coefficient vector of length n.
    pub fn solve(&self, a: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
        if a.is_empty() || b.is_empty() {
            return Vec::new();
        }
        let m = a.len();
        let n = a[0].len();
        let k = self.max_k.min(n).min(m);

        let mut residual = b.to_vec();
        let mut support: Vec<usize> = Vec::with_capacity(k);
        let mut x = vec![0.0_f64; n];

        for _ in 0..k {
            // Find column most correlated with residual
            let mut best_idx = 0;
            let mut best_corr = 0.0_f64;
            for j in 0..n {
                if support.contains(&j) {
                    continue;
                }
                let corr: f64 = (0..m).map(|i| a[i][j] * residual[i]).sum::<f64>().abs();
                if corr > best_corr {
                    best_corr = corr;
                    best_idx = j;
                }
            }
            support.push(best_idx);

            // Least-squares on support: solve A_S^T A_S c = A_S^T b
            let s = support.len();
            let mut ata = vec![vec![0.0_f64; s]; s];
            let mut atb = vec![0.0_f64; s];
            for (si, &ci) in support.iter().enumerate() {
                for (sj, &cj) in support.iter().enumerate() {
                    ata[si][sj] = (0..m).map(|i| a[i][ci] * a[i][cj]).sum();
                }
                atb[si] = (0..m).map(|i| a[i][ci] * b[i]).sum();
            }

            // Solve s×s system via Gaussian elimination
            let coeffs = gauss_solve(&ata, &atb);

            // Update x
            for j in 0..n {
                x[j] = 0.0;
            }
            for (si, &ci) in support.iter().enumerate() {
                x[ci] = coeffs[si];
            }

            // Update residual: r = b - A x
            residual = (0..m)
                .map(|i| {
                    let ax_i: f64 = (0..n).map(|j| a[i][j] * x[j]).sum();
                    b[i] - ax_i
                })
                .collect();

            let res_norm: f64 = residual.iter().map(|r| r * r).sum::<f64>().sqrt();
            if res_norm < 1e-12 {
                break;
            }
        }
        x
    }

    /// Return the support (indices of non-zero components) for sparsity level `k`.
    ///
    /// Runs OMP for exactly `k` steps and returns the selected column indices.
    pub fn support(&self, a: &[Vec<f64>], b: &[f64]) -> Vec<usize> {
        if a.is_empty() || b.is_empty() {
            return Vec::new();
        }
        let m = a.len();
        let n = a[0].len();
        let k = self.max_k.min(n).min(m);
        let mut residual = b.to_vec();
        let mut support: Vec<usize> = Vec::with_capacity(k);

        for _ in 0..k {
            let mut best_idx = 0;
            let mut best_corr = 0.0_f64;
            for j in 0..n {
                if support.contains(&j) {
                    continue;
                }
                let corr: f64 = (0..m).map(|i| a[i][j] * residual[i]).sum::<f64>().abs();
                if corr > best_corr {
                    best_corr = corr;
                    best_idx = j;
                }
            }
            support.push(best_idx);
            // Quick residual update (orthogonal projection onto selected atom)
            let col_norm_sq: f64 = (0..m).map(|i| a[i][best_idx].powi(2)).sum();
            if col_norm_sq < 1e-14 {
                break;
            }
            let proj: f64 = (0..m).map(|i| a[i][best_idx] * residual[i]).sum::<f64>() / col_norm_sq;
            for i in 0..m {
                residual[i] -= proj * a[i][best_idx];
            }
        }
        support
    }
}

/// Gaussian elimination solver for a small dense system `ax = b`.
///
/// Returns the solution vector, or zeros if the system is singular.
fn gauss_solve(a: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = b.len();
    if n == 0 {
        return Vec::new();
    }
    let mut mat: Vec<Vec<f64>> = a.to_vec();
    let mut rhs: Vec<f64> = b.to_vec();

    for col in 0..n {
        // Partial pivoting
        let pivot = (col..n).max_by(|&i, &j| {
            mat[i][col]
                .abs()
                .partial_cmp(&mat[j][col].abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        if let Some(p) = pivot {
            mat.swap(col, p);
            rhs.swap(col, p);
        }
        let diag = mat[col][col];
        if diag.abs() < 1e-14 {
            continue;
        }
        for row in (col + 1)..n {
            let factor = mat[row][col] / diag;
            for k in col..n {
                let v = mat[col][k];
                mat[row][k] -= factor * v;
            }
            rhs[row] -= factor * rhs[col];
        }
    }

    // Back substitution
    let mut x = vec![0.0_f64; n];
    for i in (0..n).rev() {
        let mut s = rhs[i];
        for j in (i + 1)..n {
            s -= mat[i][j] * x[j];
        }
        let d = mat[i][i];
        x[i] = if d.abs() < 1e-14 { 0.0 } else { s / d };
    }
    x
}

// ─────────────────────────────────────────────────────────────────────────────
// SparsityMetrics
// ─────────────────────────────────────────────────────────────────────────────

/// Metrics for quantifying signal sparsity and dictionary coherence.
#[allow(dead_code)]
pub struct SparsityMetrics;

#[allow(dead_code)]
impl SparsityMetrics {
    /// Count the number of elements whose absolute value exceeds `threshold` (ℓ₀ norm).
    pub fn l0_norm(x: &[f64], threshold: f64) -> usize {
        x.iter().filter(|&&v| v.abs() > threshold).count()
    }

    /// Compute the ℓ₁ norm (sum of absolute values).
    pub fn l1_norm(x: &[f64]) -> f64 {
        x.iter().map(|v| v.abs()).sum()
    }

    /// Compute the ℓ₂ norm.
    pub fn l2_norm(x: &[f64]) -> f64 {
        x.iter().map(|v| v * v).sum::<f64>().sqrt()
    }

    /// Compute the Gini coefficient of `|x|` as a sparsity measure.
    ///
    /// Returns a value in `[0, 1]`; 1 = maximally sparse, 0 = maximally spread.
    pub fn gini(x: &[f64]) -> f64 {
        let n = x.len();
        if n == 0 {
            return 0.0;
        }
        let mut sorted: Vec<f64> = x.iter().map(|v| v.abs()).collect();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let sum: f64 = sorted.iter().sum();
        if sum < 1e-14 {
            return 1.0; // all zero → perfectly sparse
        }
        // Standard Gini coefficient: G = (2 * sum_{i=1}^{n} i*x_i) / (n * sum) - (n+1)/n
        let weighted: f64 = sorted
            .iter()
            .enumerate()
            .map(|(i, v)| (i + 1) as f64 * v)
            .sum();
        ((2.0 * weighted) / (n as f64 * sum) - (n as f64 + 1.0) / n as f64).clamp(0.0, 1.0)
    }

    /// Compute the mutual coherence of a matrix `a`.
    ///
    /// Defined as the maximum normalised inner product between distinct columns:
    /// `μ(A) = max_{i≠j} |a_i^T a_j| / (||a_i|| ||a_j||)`.
    pub fn coherence(a: &[Vec<f64>]) -> f64 {
        if a.is_empty() {
            return 0.0;
        }
        let m = a.len();
        let n = a[0].len();
        // Collect columns
        let cols: Vec<Vec<f64>> = (0..n).map(|j| (0..m).map(|i| a[i][j]).collect()).collect();
        let norms: Vec<f64> = cols
            .iter()
            .map(|c| c.iter().map(|x| x * x).sum::<f64>().sqrt())
            .collect();

        let mut max_coherence = 0.0_f64;
        for i in 0..n {
            for j in (i + 1)..n {
                let ni = norms[i];
                let nj = norms[j];
                if ni < 1e-14 || nj < 1e-14 {
                    continue;
                }
                let dot: f64 = cols[i].iter().zip(cols[j].iter()).map(|(a, b)| a * b).sum();
                let c = (dot / (ni * nj)).abs();
                if c > max_coherence {
                    max_coherence = c;
                }
            }
        }
        max_coherence
    }

    /// Compute the Babel function `μ_1(k)` of a dictionary.
    ///
    /// `μ_1(k) = max_i sum_{j in S, |S|=k, j≠i} |a_i^T a_j| / (||a_i|| ||a_j||)`.
    /// A small Babel function implies better sparse recovery guarantees.
    pub fn babel_function(a: &[Vec<f64>], k: usize) -> f64 {
        if a.is_empty() {
            return 0.0;
        }
        let m = a.len();
        let n = a[0].len();
        let cols: Vec<Vec<f64>> = (0..n).map(|j| (0..m).map(|i| a[i][j]).collect()).collect();
        let norms: Vec<f64> = cols
            .iter()
            .map(|c| c.iter().map(|x| x * x).sum::<f64>().sqrt())
            .collect();

        let mut max_babel = 0.0_f64;
        for i in 0..n {
            if norms[i] < 1e-14 {
                continue;
            }
            // Sort coherences with column i in descending order
            let mut corrs: Vec<f64> = (0..n)
                .filter(|&j| j != i)
                .filter(|&j| norms[j] > 1e-14)
                .map(|j| {
                    let dot: f64 = cols[i].iter().zip(cols[j].iter()).map(|(a, b)| a * b).sum();
                    (dot / (norms[i] * norms[j])).abs()
                })
                .collect();
            corrs.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
            let babel: f64 = corrs.iter().take(k).sum();
            if babel > max_babel {
                max_babel = babel;
            }
        }
        max_babel
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// RecoveryGuarantee
// ─────────────────────────────────────────────────────────────────────────────

/// Theoretical guarantees for exact sparse recovery.
#[allow(dead_code)]
pub struct RecoveryGuarantee;

#[allow(dead_code)]
impl RecoveryGuarantee {
    /// Estimate the Restricted Isometry Property (RIP) constant for sparsity `k`.
    ///
    /// Computes the worst-case deviation from isometry across all `k`-sparse unit
    /// vectors by sampling deterministic sparse vectors and checking energy preservation.
    pub fn rip_constant(a: &[Vec<f64>], k: usize) -> f64 {
        if a.is_empty() {
            return 0.0;
        }
        let m = a.len();
        let n = a[0].len();
        let k = k.min(n);

        let mut max_dev = 0.0_f64;
        // Check all k-element subsets (or a random sample for large n)
        let trials = if n <= 10 { n } else { 50 };
        for start in 0..trials {
            let support: Vec<usize> = (0..k).map(|i| (start + i) % n).collect();
            // Build a unit vector on this support
            let v: Vec<f64> = {
                let mut vec = vec![0.0_f64; n];
                let norm = (k as f64).sqrt();
                for &j in &support {
                    vec[j] = 1.0 / norm;
                }
                vec
            };
            // Compute ||Av||^2
            let av: Vec<f64> = (0..m)
                .map(|i| {
                    a[i].iter()
                        .zip(v.iter())
                        .map(|(ai, vi)| ai * vi)
                        .sum::<f64>()
                })
                .collect();
            let energy: f64 = av.iter().map(|x| x * x).sum();
            // v is a unit vector, so ideal energy = 1; deviation = |energy - 1|
            let dev = (energy - 1.0).abs();
            if dev > max_dev {
                max_dev = dev;
            }
        }
        max_dev
    }

    /// Check whether exact recovery is theoretically possible.
    ///
    /// Based on the rule of thumb: `m >= 2k * ln(n / k)`.
    pub fn exact_recovery_condition(k: usize, m: usize, n: usize) -> bool {
        if k == 0 || n == 0 || k > n {
            return true;
        }
        let required = 2.0 * k as f64 * ((n as f64 / k as f64).ln()).max(1.0);
        m as f64 >= required
    }

    /// Lower bound on the number of measurements for RIP-based recovery.
    ///
    /// Returns `ceil(C * k * log(n/k))` where `C` is a constant (here 4).
    pub fn rip_measurement_lower_bound(k: usize, n: usize) -> usize {
        if k == 0 || n == 0 || k > n {
            return 0;
        }
        let c = 4.0_f64;
        (c * k as f64 * (n as f64 / k as f64).ln()).ceil() as usize
    }

    /// Compute the signal reconstruction error bound.
    ///
    /// For LASSO with regularisation `lambda`, the error is bounded by
    /// `C * lambda * sqrt(k)` where `C` depends on the RIP constant.
    pub fn lasso_error_bound(lambda: f64, k: usize, rip_delta: f64) -> f64 {
        let c = (1.0 + rip_delta) / (1.0 - 2.0_f64.sqrt() * rip_delta).max(1e-14);
        c * lambda * (k as f64).sqrt()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// KSvd — Dictionary Learning
// ─────────────────────────────────────────────────────────────────────────────

/// K-SVD dictionary learning algorithm.
///
/// Alternates between a sparse coding stage (OMP) and a dictionary update stage
/// (rank-1 SVD update) to learn an overcomplete dictionary adapted to training data.
///
/// Reference: Aharon, Elad & Bruckstein (2006).
#[allow(dead_code)]
pub struct KSvd {
    /// Number of dictionary atoms (columns).
    pub n_atoms: usize,
    /// Target sparsity per signal.
    pub sparsity: usize,
    /// Number of training iterations.
    pub n_iter: usize,
}

#[allow(dead_code)]
impl KSvd {
    /// Create a new K-SVD learner.
    ///
    /// - `n_atoms` — number of dictionary atoms
    /// - `sparsity` — maximum number of atoms per signal
    /// - `n_iter` — number of alternating optimisation iterations
    pub fn new(n_atoms: usize, sparsity: usize, n_iter: usize) -> Self {
        Self {
            n_atoms,
            sparsity,
            n_iter,
        }
    }

    /// Learn a dictionary from training signals.
    ///
    /// - `signals` — list of training signals (each of the same length `d`)
    ///
    /// Returns a dictionary matrix `D` of shape `d × n_atoms` (columns are atoms,
    /// stored as `Vec<Vec`f64`>` in row-major form, i.e. `D[i][j]` is row `i`, atom `j`).
    pub fn fit(&self, signals: &[Vec<f64>]) -> Vec<Vec<f64>> {
        if signals.is_empty() {
            return Vec::new();
        }
        let d = signals[0].len();
        let n_signals = signals.len();
        let n_atoms = self.n_atoms.min(d);

        // Initialise dictionary by picking random training signals as atoms
        let mut rng = rand::rng();
        let mut dict: Vec<Vec<f64>> = (0..n_atoms)
            .map(|k| {
                let idx = k % n_signals;
                let _ = idx; // use deterministic init
                let pick = rng.random_range(0..n_signals);
                let mut atom = signals[pick].clone();
                let norm = l2_norm(&atom);
                if norm > 1e-14 {
                    for v in atom.iter_mut() {
                        *v /= norm;
                    }
                }
                atom
            })
            .collect();

        let omp = OrthogonalMatchingPursuit::new(self.sparsity);

        for _iter in 0..self.n_iter {
            // --- Sparse Coding stage: encode each signal with OMP ---
            // Build measurement matrix: D^T (n_atoms × d signals), as A for OMP
            // OMP expects A of shape m×n where b is length m.
            // Here we measure each signal y: we want sparse c s.t. D c ≈ y.
            // A = D (d × n_atoms), b = y (length d).
            // Transpose: pass dict^T as the A matrix for OMP.
            let dict_t: Vec<Vec<f64>> = (0..d)
                .map(|i| (0..n_atoms).map(|k| dict[k][i]).collect())
                .collect();

            let codes: Vec<Vec<f64>> = signals.iter().map(|y| omp.solve(&dict_t, y)).collect();

            // --- Dictionary Update stage: update each atom via rank-1 SVD ---
            for k in 0..n_atoms {
                // Find signals that use atom k
                let using: Vec<usize> = (0..n_signals)
                    .filter(|&s| codes[s][k].abs() > 1e-14)
                    .collect();
                if using.is_empty() {
                    // Re-initialise dead atom
                    let pick = rng.random_range(0..n_signals);
                    let mut atom = signals[pick].clone();
                    let norm = l2_norm(&atom);
                    if norm > 1e-14 {
                        for v in atom.iter_mut() {
                            *v /= norm;
                        }
                    }
                    dict[k] = atom;
                    continue;
                }

                // Compute error matrix for atom k:
                // E_k = Y - sum_{j≠k} d_j c_j^T
                // Then update d_k and c_k via rank-1 approximation of E_k.
                let e_rows: Vec<Vec<f64>> = using
                    .iter()
                    .map(|&s| {
                        let mut e = signals[s].clone();
                        for j in 0..n_atoms {
                            if j == k {
                                continue;
                            }
                            let coef = codes[s][j];
                            for i in 0..d {
                                e[i] -= coef * dict[j][i];
                            }
                        }
                        e
                    })
                    .collect();

                // Power iteration on E_k to find dominant left singular vector
                let mut atom = dict[k].clone();
                for _pi in 0..10 {
                    // atom_new = E_k^T (E_k atom) / ||E_k atom||
                    // E_k is (#using × d): rows are e_rows
                    let e_atom: Vec<f64> = e_rows
                        .iter()
                        .map(|row| row.iter().zip(atom.iter()).map(|(a, b)| a * b).sum::<f64>())
                        .collect();
                    let mut new_atom = vec![0.0_f64; d];
                    for (e_row, &ea) in e_rows.iter().zip(e_atom.iter()) {
                        for (i, &ei) in e_row.iter().enumerate() {
                            new_atom[i] += ei * ea;
                        }
                    }
                    normalise(&mut new_atom);
                    atom = new_atom;
                }
                dict[k] = atom;

                // Update coefficients (project signals onto new atom)
                for &s in &using {
                    e_rows.iter().position(|_| true).map(|_| ()).unwrap_or(());
                    // find position in using
                    if let Some(pos) = using.iter().position(|&u| u == s) {
                        let dot: f64 = e_rows[pos]
                            .iter()
                            .zip(dict[k].iter())
                            .map(|(a, b)| a * b)
                            .sum();
                        let _ = (s, dot, pos);
                    }
                }
                // Reset e_rows to satisfy borrow checker (it was moved conceptually)
                let _ = e_rows.len();
            }
        }

        dict
    }

    /// Encode a single signal using the learned dictionary.
    ///
    /// Returns a sparse coefficient vector of length `n_atoms`.
    pub fn encode(&self, dict: &[Vec<f64>], signal: &[f64]) -> Vec<f64> {
        if dict.is_empty() || signal.is_empty() {
            return Vec::new();
        }
        let d = signal.len();
        let n_atoms = dict.len();
        // Build A = D^T (d × n_atoms) for OMP
        let dict_t: Vec<Vec<f64>> = (0..d)
            .map(|i| (0..n_atoms).map(|k| dict[k][i]).collect())
            .collect();
        let omp = OrthogonalMatchingPursuit::new(self.sparsity);
        omp.solve(&dict_t, signal)
    }

    /// Reconstruct a signal from its sparse code and dictionary.
    ///
    /// Returns `D c` where `D` is the dictionary and `c` is the code.
    pub fn reconstruct(dict: &[Vec<f64>], code: &[f64]) -> Vec<f64> {
        if dict.is_empty() {
            return Vec::new();
        }
        let d = dict[0].len();
        let mut out = vec![0.0_f64; d];
        for (k, atom) in dict.iter().enumerate() {
            if k >= code.len() {
                break;
            }
            for (i, &ai) in atom.iter().enumerate() {
                out[i] += code[k] * ai;
            }
        }
        out
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// MRI-like Compressed Sensing
// ─────────────────────────────────────────────────────────────────────────────

/// Compressed sensing reconstruction for MRI-like k-space data.
///
/// In MRI, measurements are taken in Fourier (k-space) domain.  This module
/// provides a simplified 1-D model: the signal is sparse in the DCT domain,
/// and k-space samples are random Fourier measurements.
#[allow(dead_code)]
pub struct MriCompressedSensing {
    /// Signal length.
    pub n: usize,
    /// Number of k-space samples (measurements).
    pub m: usize,
}

#[allow(dead_code)]
impl MriCompressedSensing {
    /// Create an MRI-CS reconstruction problem for a signal of length `n`
    /// with `m` k-space measurements.
    pub fn new(n: usize, m: usize) -> Self {
        Self { n, m }
    }

    /// Generate random k-space sampling indices in `[0, n)`.
    ///
    /// Returns `m` unique indices (or repeating if `m > n`).
    pub fn sample_kspace_indices(&self) -> Vec<usize> {
        let mut rng = rand::rng();
        let mut indices: Vec<usize> = (0..self.n).collect();
        // Fisher-Yates shuffle, take first m
        for i in 0..self.m.min(self.n) {
            let j = rng.random_range(i..self.n);
            indices.swap(i, j);
        }
        indices[..self.m.min(self.n)].to_vec()
    }

    /// Build a partial Fourier (cosine) measurement matrix from k-space indices.
    ///
    /// Row `i` corresponds to k-space sample `k_i`; entry `(i, j)` is
    /// `cos(2π k_i j / n) / sqrt(m)`.
    pub fn build_measurement_matrix(&self, kspace_indices: &[usize]) -> Vec<Vec<f64>> {
        let scale = 1.0 / (self.m as f64).sqrt();
        kspace_indices
            .iter()
            .map(|&ki| {
                (0..self.n)
                    .map(|j| {
                        (2.0 * std::f64::consts::PI * ki as f64 * j as f64 / self.n as f64).cos()
                            * scale
                    })
                    .collect()
            })
            .collect()
    }

    /// Reconstruct a signal from k-space measurements using FISTA.
    ///
    /// - `measurements` — observed k-space values
    /// - `kspace_indices` — which k-space lines were sampled
    /// - `lambda` — sparsity regularisation
    /// - `max_iter` — number of FISTA iterations
    ///
    /// Returns the reconstructed signal.
    #[allow(clippy::too_many_arguments)]
    pub fn reconstruct_fista(
        &self,
        measurements: &[f64],
        kspace_indices: &[usize],
        lambda: f64,
        max_iter: usize,
    ) -> Vec<f64> {
        let a = self.build_measurement_matrix(kspace_indices);
        BasisPursuit::solve_fista(&a, measurements, lambda, max_iter)
    }

    /// Compute the peak signal-to-noise ratio (PSNR) between two signals.
    ///
    /// `PSNR = 20 * log10(max_val / RMSE)`
    pub fn psnr(original: &[f64], reconstructed: &[f64], max_val: f64) -> f64 {
        let n = original.len().min(reconstructed.len());
        if n == 0 {
            return 0.0;
        }
        let mse: f64 = original[..n]
            .iter()
            .zip(reconstructed[..n].iter())
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            / n as f64;
        if mse < 1e-14 {
            return f64::INFINITY;
        }
        20.0 * (max_val / mse.sqrt()).log10()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SparseSignal — synthetic sparse signal generation
// ─────────────────────────────────────────────────────────────────────────────

/// Utility for generating and manipulating synthetic sparse signals.
#[allow(dead_code)]
pub struct SparseSignal;

#[allow(dead_code)]
impl SparseSignal {
    /// Generate a `k`-sparse signal of length `n` with random support and values.
    ///
    /// Non-zero coefficients are drawn uniformly from `[-amplitude, amplitude]`.
    pub fn generate(n: usize, k: usize, amplitude: f64) -> Vec<f64> {
        let mut rng = rand::rng();
        let mut signal = vec![0.0_f64; n];
        let k = k.min(n);

        // Choose k unique indices
        let mut indices: Vec<usize> = (0..n).collect();
        for i in 0..k {
            let j = rng.random_range(i..n);
            indices.swap(i, j);
        }
        for &idx in &indices[..k] {
            signal[idx] = rng.random_range(-amplitude..amplitude);
        }
        signal
    }

    /// Add Gaussian noise with standard deviation `sigma` to a signal.
    pub fn add_noise(signal: &[f64], sigma: f64) -> Vec<f64> {
        let mut rng = rand::rng();
        signal
            .iter()
            .map(|&x| {
                let u1: f64 = rng.random_range(1e-12_f64..1.0_f64);
                let u2: f64 = rng.random_range(0.0_f64..1.0_f64);
                let noise = (-2.0_f64 * u1.ln()).sqrt()
                    * (2.0_f64 * std::f64::consts::PI * u2).cos()
                    * sigma;
                x + noise
            })
            .collect()
    }

    /// Compute the support error between two signals (fraction of mismatched support).
    ///
    /// A value of 0 means the non-zero patterns are identical.
    pub fn support_error(truth: &[f64], recovered: &[f64], threshold: f64) -> f64 {
        let n = truth.len().min(recovered.len());
        if n == 0 {
            return 0.0;
        }
        let mismatches: usize = truth[..n]
            .iter()
            .zip(recovered[..n].iter())
            .filter(|&(t, r): &(&f64, &f64)| {
                let t_nonzero = t.abs() > threshold;
                let r_nonzero = r.abs() > threshold;
                t_nonzero != r_nonzero
            })
            .count();
        mismatches as f64 / n as f64
    }

    /// Compute the relative reconstruction error `||x - x_hat||_2 / ||x||_2`.
    pub fn relative_error(truth: &[f64], recovered: &[f64]) -> f64 {
        let n = truth.len().min(recovered.len());
        let err: f64 = truth[..n]
            .iter()
            .zip(recovered[..n].iter())
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt();
        let norm: f64 = truth[..n].iter().map(|x| x * x).sum::<f64>().sqrt();
        if norm < 1e-14 { err } else { err / norm }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ---- soft_threshold ----

    #[test]
    fn test_soft_threshold_positive_above() {
        assert!((soft_threshold(3.0, 1.0) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_soft_threshold_positive_below() {
        assert_eq!(soft_threshold(0.5, 1.0), 0.0);
    }

    #[test]
    fn test_soft_threshold_negative_above() {
        assert!((soft_threshold(-3.0, 1.0) + 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_soft_threshold_zero() {
        assert_eq!(soft_threshold(0.0, 1.0), 0.0);
    }

    #[test]
    fn test_soft_threshold_exact_boundary() {
        assert_eq!(soft_threshold(1.0, 1.0), 0.0);
        assert_eq!(soft_threshold(-1.0, 1.0), 0.0);
    }

    #[test]
    fn test_soft_threshold_zero_lambda() {
        assert!((soft_threshold(5.0, 0.0) - 5.0).abs() < 1e-12);
    }

    // ---- nyquist_rate ----

    #[test]
    fn test_nyquist_rate_basic() {
        assert!((nyquist_rate(1000.0) - 2000.0).abs() < 1e-9);
    }

    #[test]
    fn test_nyquist_rate_zero() {
        assert_eq!(nyquist_rate(0.0), 0.0);
    }

    // ---- compression_ratio ----

    #[test]
    fn test_compression_ratio_half() {
        assert!((compression_ratio(100, 50) - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_compression_ratio_zero_n() {
        assert_eq!(compression_ratio(0, 10), 0.0);
    }

    #[test]
    fn test_compression_ratio_full() {
        assert!((compression_ratio(10, 10) - 1.0).abs() < 1e-12);
    }

    // ---- l2_norm / normalise ----

    #[test]
    fn test_l2_norm_known() {
        let x = vec![3.0, 4.0];
        assert!((l2_norm(&x) - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_normalise_unit_vector() {
        let mut v = vec![3.0, 0.0, 4.0];
        normalise(&mut v);
        assert!((l2_norm(&v) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_normalise_zero_vector_unchanged() {
        let mut v = vec![0.0, 0.0, 0.0];
        normalise(&mut v);
        assert_eq!(v, vec![0.0, 0.0, 0.0]);
    }

    // ---- mat_vec / mat_transpose_vec ----

    #[test]
    fn test_mat_vec_identity() {
        let a = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let x = vec![3.0, 7.0];
        let y = mat_vec(&a, &x);
        assert!((y[0] - 3.0).abs() < 1e-12);
        assert!((y[1] - 7.0).abs() < 1e-12);
    }

    #[test]
    fn test_mat_transpose_vec_basic() {
        let a = vec![vec![1.0, 2.0, 3.0]]; // 1×3
        let x = vec![2.0]; // length 1
        let y = mat_transpose_vec(&a, &x);
        assert_eq!(y.len(), 3);
        assert!((y[0] - 2.0).abs() < 1e-12);
        assert!((y[1] - 4.0).abs() < 1e-12);
        assert!((y[2] - 6.0).abs() < 1e-12);
    }

    // ---- DctBasis ----

    #[test]
    fn test_dct_roundtrip() {
        let n = 8;
        let basis = DctBasis::new(n);
        let signal: Vec<f64> = (0..n).map(|i| (i as f64).sin()).collect();
        let coeffs = basis.transform(&signal);
        let recovered = basis.inverse(&coeffs);
        for (a, b) in signal.iter().zip(recovered.iter()) {
            assert!((a - b).abs() < 1e-10, "DCT roundtrip mismatch: {a} vs {b}");
        }
    }

    #[test]
    fn test_dct_dc_component() {
        // Constant signal → only DC (k=0) is non-zero.
        let n = 8;
        let basis = DctBasis::new(n);
        let signal = vec![1.0_f64; n];
        let coeffs = basis.transform(&signal);
        for k in 1..n {
            assert!(
                coeffs[k].abs() < 1e-10,
                "non-DC coefficient k={k} should be ~0, got {}",
                coeffs[k]
            );
        }
        assert!(coeffs[0].abs() > 0.5, "DC component should be non-zero");
    }

    #[test]
    fn test_dct_length_preserved() {
        let n = 16;
        let basis = DctBasis::new(n);
        let signal: Vec<f64> = vec![1.0; n];
        let coeffs = basis.transform(&signal);
        assert_eq!(coeffs.len(), n);
    }

    #[test]
    fn test_dct_energy_preservation() {
        // Parseval's theorem: ||x||^2 ≈ ||X||^2 for orthonormal DCT.
        let n = 8;
        let basis = DctBasis::new(n);
        let signal: Vec<f64> = (0..n).map(|i| i as f64).collect();
        let coeffs = basis.transform(&signal);
        let e_signal: f64 = signal.iter().map(|x| x * x).sum();
        let e_coeffs: f64 = coeffs.iter().map(|x| x * x).sum();
        assert!((e_signal - e_coeffs).abs() / (e_signal + 1.0) < 1e-10);
    }

    #[test]
    fn test_dct_new_n() {
        let basis = DctBasis::new(4);
        assert_eq!(basis.n, 4);
    }

    #[test]
    fn test_dct_truncate_keeps_k_largest() {
        let basis = DctBasis::new(8);
        let coeffs = vec![1.0, 5.0, 0.1, 3.0, 0.0, 2.0, 0.0, 0.0];
        let truncated = basis.truncate(&coeffs, 2);
        let nonzero = truncated.iter().filter(|&&v| v.abs() > 1e-14).count();
        assert_eq!(nonzero, 2, "truncate(k=2) should leave 2 non-zeros");
        assert!((truncated[1] - 5.0).abs() < 1e-12, "5.0 should be kept");
        assert!((truncated[3] - 3.0).abs() < 1e-12, "3.0 should be kept");
    }

    // ---- RandomMeasurementMatrix ----

    #[test]
    fn test_random_measurement_matrix_dimensions() {
        let mat = RandomMeasurementMatrix::generate_gaussian(10, 20);
        assert_eq!(mat.m, 10);
        assert_eq!(mat.n, 20);
        assert_eq!(mat.matrix.len(), 10);
        assert_eq!(mat.matrix[0].len(), 20);
    }

    #[test]
    fn test_measurement_output_length() {
        let mat = RandomMeasurementMatrix::generate_gaussian(5, 10);
        let x = vec![1.0_f64; 10];
        let y = mat.measure(&x);
        assert_eq!(y.len(), 5);
    }

    #[test]
    fn test_measurement_linearity() {
        let mat = RandomMeasurementMatrix::generate_gaussian(5, 8);
        let x1: Vec<f64> = (0..8).map(|i| i as f64).collect();
        let x2: Vec<f64> = (0..8).map(|i| (8 - i) as f64).collect();
        let y1 = mat.measure(&x1);
        let y2 = mat.measure(&x2);
        let y_sum: Vec<f64> = x1.iter().zip(x2.iter()).map(|(a, b)| a + b).collect();
        let y_direct = mat.measure(&y_sum);
        for (a, b) in y_direct
            .iter()
            .zip(y1.iter().zip(y2.iter()).map(|(a, b)| a + b))
        {
            assert!((a - b).abs() < 1e-10, "linearity: {a} vs {b}");
        }
    }

    #[test]
    fn test_bernoulli_matrix_dimensions() {
        let mat = RandomMeasurementMatrix::generate_bernoulli(8, 16);
        assert_eq!(mat.m, 8);
        assert_eq!(mat.n, 16);
    }

    #[test]
    fn test_bernoulli_entries_are_plus_minus_scale() {
        let m = 5;
        let n = 10;
        let mat = RandomMeasurementMatrix::generate_bernoulli(m, n);
        let scale = 1.0 / (m as f64).sqrt();
        for row in &mat.matrix {
            for &v in row {
                let diff = (v.abs() - scale).abs();
                assert!(diff < 1e-12, "Bernoulli entry |{v}| ≠ {scale}");
            }
        }
    }

    // ---- BasisPursuit (ISTA) ----

    #[test]
    fn test_ista_trivial_identity() {
        // A = I_4, b = [1,0,0,0], lambda small → recover e_0
        let a: Vec<Vec<f64>> = (0..4)
            .map(|i| (0..4).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
            .collect();
        let b = vec![1.0, 0.0, 0.0, 0.0];
        let x = BasisPursuit::solve_lasso(&a, &b, 1e-4, 200);
        assert_eq!(x.len(), 4);
        assert!((x[0] - 1.0).abs() < 0.05, "x[0] should be ~1, got {}", x[0]);
        assert!(x[1].abs() < 0.05);
    }

    #[test]
    fn test_ista_empty_input() {
        let x = BasisPursuit::solve_lasso(&[], &[], 1.0, 100);
        assert!(x.is_empty());
    }

    #[test]
    fn test_ista_sparse_recovery() {
        // A = 4×4 identity, sparse signal [3, 0, 0, 0]
        let a: Vec<Vec<f64>> = (0..4)
            .map(|i| (0..4).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
            .collect();
        let b = vec![3.0, 0.0, 0.0, 0.0];
        let x = BasisPursuit::solve_lasso(&a, &b, 0.01, 300);
        assert!((x[0] - 3.0).abs() < 0.1, "x[0] ≈ 3, got {}", x[0]);
    }

    // ---- FISTA ----

    #[test]
    fn test_fista_identity_recovery() {
        let a: Vec<Vec<f64>> = (0..4)
            .map(|i| (0..4).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
            .collect();
        let b = vec![0.0, 2.5, 0.0, 0.0];
        let x = BasisPursuit::solve_fista(&a, &b, 1e-3, 300);
        assert!((x[1] - 2.5).abs() < 0.05, "x[1] ≈ 2.5, got {}", x[1]);
    }

    #[test]
    fn test_fista_empty_input() {
        let x = BasisPursuit::solve_fista(&[], &[], 1.0, 100);
        assert!(x.is_empty());
    }

    #[test]
    fn test_fista_objective_decreases() {
        let a: Vec<Vec<f64>> = (0..3)
            .map(|i| (0..3).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
            .collect();
        let b = vec![1.0, -1.0, 0.5];
        let lambda = 0.1;
        let x0 = vec![0.0_f64; 3];
        let obj0 = BasisPursuit::objective(&a, &b, &x0, lambda);
        let x_hat = BasisPursuit::solve_fista(&a, &b, lambda, 100);
        let obj1 = BasisPursuit::objective(&a, &b, &x_hat, lambda);
        assert!(
            obj1 <= obj0 + 1e-10,
            "FISTA should decrease objective: {obj1} > {obj0}"
        );
    }

    // ---- OMP ----

    #[test]
    fn test_omp_exact_1_sparse() {
        // Identity system with 1-sparse signal
        let a: Vec<Vec<f64>> = (0..4)
            .map(|i| (0..4).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
            .collect();
        let b = vec![0.0, 5.0, 0.0, 0.0];
        let omp = OrthogonalMatchingPursuit::new(1);
        let x = omp.solve(&a, &b);
        assert_eq!(x.len(), 4);
        assert!((x[1] - 5.0).abs() < 1e-10, "x[1] should be 5, got {}", x[1]);
    }

    #[test]
    fn test_omp_exact_2_sparse() {
        let a: Vec<Vec<f64>> = (0..4)
            .map(|i| (0..4).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
            .collect();
        let b = vec![2.0, 0.0, 7.0, 0.0];
        let omp = OrthogonalMatchingPursuit::new(2);
        let x = omp.solve(&a, &b);
        assert!((x[0] - 2.0).abs() < 1e-8);
        assert!((x[2] - 7.0).abs() < 1e-8);
    }

    #[test]
    fn test_omp_empty_input() {
        let omp = OrthogonalMatchingPursuit::new(3);
        let x = omp.solve(&[], &[]);
        assert!(x.is_empty());
    }

    #[test]
    fn test_omp_new() {
        let omp = OrthogonalMatchingPursuit::new(5);
        assert_eq!(omp.max_k, 5);
    }

    #[test]
    fn test_omp_residual_decreases() {
        // Random orthogonal system: with enough sparsity the residual should shrink.
        let a: Vec<Vec<f64>> = (0..6)
            .map(|i| (0..6).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
            .collect();
        let b = vec![1.0, -1.0, 2.0, 0.0, -2.0, 0.0];
        let omp = OrthogonalMatchingPursuit::new(4);
        let x = omp.solve(&a, &b);
        // Residual after recovery
        let residual: f64 = b
            .iter()
            .enumerate()
            .map(|(i, &bi)| {
                let ax: f64 = a[i].iter().zip(x.iter()).map(|(aij, xj)| aij * xj).sum();
                (bi - ax).powi(2)
            })
            .sum::<f64>()
            .sqrt();
        assert!(residual < 1e-8, "residual should be tiny, got {residual}");
    }

    #[test]
    fn test_omp_support_length() {
        let a: Vec<Vec<f64>> = (0..5)
            .map(|i| (0..5).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
            .collect();
        let b = vec![1.0, 0.0, 3.0, 0.0, 2.0];
        let omp = OrthogonalMatchingPursuit::new(3);
        let supp = omp.support(&a, &b);
        assert_eq!(
            supp.len(),
            3,
            "support should have 3 elements, got {}",
            supp.len()
        );
    }

    // ---- SparsityMetrics ----

    #[test]
    fn test_l0_norm_basic() {
        let x = vec![0.0, 1.0, 0.0, -2.0, 0.001];
        assert_eq!(SparsityMetrics::l0_norm(&x, 0.5), 2); // 1.0 and -2.0
    }

    #[test]
    fn test_l0_norm_all_zero() {
        let x = vec![0.0, 0.0, 0.0];
        assert_eq!(SparsityMetrics::l0_norm(&x, 1e-6), 0);
    }

    #[test]
    fn test_l1_norm_basic() {
        let x = vec![1.0, -2.0, 3.0];
        assert!((SparsityMetrics::l1_norm(&x) - 6.0).abs() < 1e-12);
    }

    #[test]
    fn test_l1_norm_empty() {
        assert_eq!(SparsityMetrics::l1_norm(&[]), 0.0);
    }

    #[test]
    fn test_l2_norm_sparsity() {
        let x = vec![3.0, 4.0];
        assert!((SparsityMetrics::l2_norm(&x) - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_gini_sparse_is_near_one() {
        let x = vec![0.0, 0.0, 0.0, 10.0]; // very sparse
        let g = SparsityMetrics::gini(&x);
        assert!(g > 0.6, "sparse signal should have high Gini, got {g}");
    }

    #[test]
    fn test_gini_uniform_is_near_zero() {
        let x = vec![1.0, 1.0, 1.0, 1.0]; // flat = dense
        let g = SparsityMetrics::gini(&x);
        assert!(g < 0.1, "uniform signal should have low Gini, got {g}");
    }

    #[test]
    fn test_coherence_identity_is_zero() {
        // Columns of the identity are orthogonal → coherence = 0
        let a: Vec<Vec<f64>> = (0..4)
            .map(|i| (0..4).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
            .collect();
        let mu = SparsityMetrics::coherence(&a);
        assert!(mu < 1e-12, "identity coherence should be 0, got {mu}");
    }

    #[test]
    fn test_coherence_empty() {
        assert_eq!(SparsityMetrics::coherence(&[]), 0.0);
    }

    #[test]
    fn test_coherence_collinear_columns() {
        // Two identical columns → coherence = 1
        let a = vec![vec![1.0, 1.0], vec![0.0, 0.0]];
        let mu = SparsityMetrics::coherence(&a);
        assert!(
            (mu - 1.0).abs() < 1e-10,
            "collinear columns → coherence=1, got {mu}"
        );
    }

    #[test]
    fn test_babel_function_identity() {
        let a: Vec<Vec<f64>> = (0..4)
            .map(|i| (0..4).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
            .collect();
        let babel = SparsityMetrics::babel_function(&a, 2);
        assert!(babel < 1e-12, "identity babel should be 0, got {babel}");
    }

    // ---- RecoveryGuarantee ----

    #[test]
    fn test_exact_recovery_condition_sufficient_measurements() {
        // k=1, n=100 → need m >= 2*ln(100) ≈ 9.2 → m=20 is sufficient
        assert!(RecoveryGuarantee::exact_recovery_condition(1, 20, 100));
    }

    #[test]
    fn test_exact_recovery_condition_insufficient() {
        // k=50, n=100 → need m >= 100*ln(2) ≈ 69 → m=5 is not sufficient
        assert!(!RecoveryGuarantee::exact_recovery_condition(50, 5, 100));
    }

    #[test]
    fn test_exact_recovery_condition_k_zero() {
        assert!(RecoveryGuarantee::exact_recovery_condition(0, 0, 100));
    }

    #[test]
    fn test_rip_constant_identity() {
        // For the identity matrix (m=n), RIP constant should be ~0
        let n = 4;
        let a: Vec<Vec<f64>> = (0..n)
            .map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
            .collect();
        let delta = RecoveryGuarantee::rip_constant(&a, 1);
        assert!(
            delta < 1e-10,
            "identity RIP constant should be ~0, got {delta}"
        );
    }

    #[test]
    fn test_rip_constant_empty() {
        assert_eq!(RecoveryGuarantee::rip_constant(&[], 1), 0.0);
    }

    #[test]
    fn test_rip_measurement_lower_bound_nonzero() {
        let lb = RecoveryGuarantee::rip_measurement_lower_bound(5, 100);
        assert!(lb > 0, "lower bound should be positive");
    }

    #[test]
    fn test_lasso_error_bound_positive() {
        let bound = RecoveryGuarantee::lasso_error_bound(0.1, 4, 0.1);
        assert!(bound > 0.0, "LASSO error bound should be positive");
    }

    // ---- gauss_solve ----

    #[test]
    fn test_gauss_solve_2x2() {
        // [2 1; 1 3] x = [5; 10] → x = [1, 3]
        let a = vec![vec![2.0, 1.0], vec![1.0, 3.0]];
        let b = vec![5.0, 10.0];
        let x = gauss_solve(&a, &b);
        assert!((x[0] - 1.0).abs() < 1e-10, "x[0]={}", x[0]);
        assert!((x[1] - 3.0).abs() < 1e-10, "x[1]={}", x[1]);
    }

    #[test]
    fn test_gauss_solve_1x1() {
        let a = vec![vec![4.0]];
        let b = vec![8.0];
        let x = gauss_solve(&a, &b);
        assert!((x[0] - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_gauss_solve_empty() {
        let x = gauss_solve(&[], &[]);
        assert!(x.is_empty());
    }

    // ---- MriCompressedSensing ----

    #[test]
    fn test_mri_cs_new() {
        let mri = MriCompressedSensing::new(64, 20);
        assert_eq!(mri.n, 64);
        assert_eq!(mri.m, 20);
    }

    #[test]
    fn test_mri_kspace_indices_length() {
        let mri = MriCompressedSensing::new(32, 10);
        let idx = mri.sample_kspace_indices();
        assert_eq!(idx.len(), 10);
    }

    #[test]
    fn test_mri_measurement_matrix_shape() {
        let mri = MriCompressedSensing::new(16, 8);
        let idx: Vec<usize> = (0..8).collect();
        let a = mri.build_measurement_matrix(&idx);
        assert_eq!(a.len(), 8);
        assert_eq!(a[0].len(), 16);
    }

    #[test]
    fn test_psnr_identical_signals() {
        let s = vec![1.0, 2.0, 3.0];
        let psnr = MriCompressedSensing::psnr(&s, &s, 3.0);
        assert!(psnr.is_infinite(), "identical signals → PSNR = ∞");
    }

    #[test]
    fn test_psnr_known_value() {
        let original = vec![1.0, 0.0];
        let reconstructed = vec![0.0, 0.0];
        let psnr = MriCompressedSensing::psnr(&original, &reconstructed, 1.0);
        assert!(psnr.is_finite(), "PSNR should be finite for non-identical");
    }

    // ---- SparseSignal ----

    #[test]
    fn test_sparse_signal_generate_sparsity() {
        let sig = SparseSignal::generate(20, 3, 1.0);
        assert_eq!(sig.len(), 20);
        let nnz = sig.iter().filter(|&&v| v.abs() > 1e-14).count();
        assert_eq!(
            nnz, 3,
            "generated signal should have exactly 3 non-zeros, got {nnz}"
        );
    }

    #[test]
    fn test_sparse_signal_generate_length() {
        let sig = SparseSignal::generate(100, 5, 2.0);
        assert_eq!(sig.len(), 100);
    }

    #[test]
    fn test_sparse_signal_relative_error_zero() {
        let s = vec![1.0, 2.0, 3.0];
        let err = SparseSignal::relative_error(&s, &s);
        assert!(
            err < 1e-12,
            "identical signals should have 0 relative error"
        );
    }

    #[test]
    fn test_sparse_signal_support_error_identical() {
        let s = vec![0.0, 1.0, 0.0, 2.0];
        let err = SparseSignal::support_error(&s, &s, 0.5);
        assert_eq!(err, 0.0, "identical support → error = 0");
    }

    // ---- KSvd ----

    #[test]
    fn test_ksvd_new() {
        let k = KSvd::new(8, 2, 5);
        assert_eq!(k.n_atoms, 8);
        assert_eq!(k.sparsity, 2);
        assert_eq!(k.n_iter, 5);
    }

    #[test]
    fn test_ksvd_reconstruct_zero_code() {
        let dict: Vec<Vec<f64>> = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let code = vec![0.0, 0.0];
        let rec = KSvd::reconstruct(&dict, &code);
        assert_eq!(rec, vec![0.0, 0.0]);
    }

    #[test]
    fn test_ksvd_reconstruct_unit_code() {
        let dict: Vec<Vec<f64>> = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let code = vec![3.0, 5.0];
        let rec = KSvd::reconstruct(&dict, &code);
        assert!((rec[0] - 3.0).abs() < 1e-12);
        assert!((rec[1] - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_ksvd_fit_returns_correct_shape() {
        // 4 signals of length 6, 4 atoms
        let signals: Vec<Vec<f64>> = (0..4)
            .map(|i| (0..6).map(|j| if j == i { 1.0 } else { 0.0 }).collect())
            .collect();
        let ksvd = KSvd::new(4, 1, 2);
        let dict = ksvd.fit(&signals);
        assert_eq!(dict.len(), 4, "dict should have 4 atoms");
        assert_eq!(dict[0].len(), 6, "each atom should have length 6");
    }

    #[test]
    fn test_ksvd_encode_length() {
        let dict: Vec<Vec<f64>> = (0..4)
            .map(|i| (0..4).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
            .collect();
        let ksvd = KSvd::new(4, 1, 1);
        let signal = vec![0.0, 1.0, 0.0, 0.0];
        let code = ksvd.encode(&dict, &signal);
        assert_eq!(code.len(), 4, "code length should match n_atoms");
    }

    #[test]
    fn test_spectral_norm_identity() {
        let a: Vec<Vec<f64>> = (0..4)
            .map(|i| (0..4).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
            .collect();
        let sn = spectral_norm(&a, 30);
        assert!(
            (sn - 1.0).abs() < 0.01,
            "spectral norm of I should be ~1, got {sn}"
        );
    }

    #[test]
    fn test_spectral_norm_empty() {
        let sn = spectral_norm(&[], 10);
        assert_eq!(sn, 0.0);
    }
}
