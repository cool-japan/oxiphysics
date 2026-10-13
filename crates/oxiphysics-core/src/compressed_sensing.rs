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

use rand::rngs::StdRng;
use rand::{RngExt, SeedableRng};

// ─────────────────────────────────────────────────────────────────────────────
// Free functions
// ─────────────────────────────────────────────────────────────────────────────

/// Apply element-wise soft-thresholding: `sign(x) * max(|x| - lambda, 0)`.
///
/// Used as a proximal operator in iterative shrinkage-thresholding algorithms.
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
pub fn nyquist_rate(bandwidth: f64) -> f64 {
    2.0 * bandwidth
}

/// Compute the compression ratio `m / n`.
///
/// A ratio less than 1 indicates sub-Nyquist sampling.
pub fn compression_ratio(n: usize, m: usize) -> f64 {
    if n == 0 {
        return 0.0;
    }
    m as f64 / n as f64
}

/// Compute the ℓ₂ norm of a slice.
///
/// Returns `sqrt(sum of squares)`.
pub fn l2_norm(x: &[f64]) -> f64 {
    x.iter().map(|v| v * v).sum::<f64>().sqrt()
}

/// Normalise a vector to unit ℓ₂ norm in place.
///
/// If the norm is smaller than `1e-14` the vector is left unchanged.
pub fn normalise(x: &mut [f64]) {
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
pub fn mat_vec(a: &[Vec<f64>], x: &[f64]) -> Vec<f64> {
    a.iter()
        .map(|row| row.iter().zip(x.iter()).map(|(ai, xi)| ai * xi).sum())
        .collect()
}

/// Compute the transposed matrix-vector product `y = A^T x`.
///
/// `a` is row-major with shape `m × n`; returns a vector of length `n`.
pub fn mat_transpose_vec(a: &[Vec<f64>], x: &[f64]) -> Vec<f64> {
    if a.is_empty() {
        return Vec::new();
    }
    let n = a[0].len();
    let mut y = vec![0.0_f64; n];
    for (row, &xi) in a.iter().zip(x.iter()) {
        for (yj, &aij) in y.iter_mut().zip(row.iter()) {
            *yj += aij * xi;
        }
    }
    y
}

/// Estimate the spectral norm (largest singular value) of `a` via power iteration.
///
/// Runs `max_iter` iterations; returns an approximation to `||A||_2`.
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
pub struct DctBasis {
    /// Length of the signal (number of samples).
    pub n: usize,
}

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
        let pi_over_n = std::f64::consts::PI / n as f64;
        (0..n)
            .map(|k| {
                let sum: f64 = x[..n]
                    .iter()
                    .enumerate()
                    .map(|(j, &xj)| xj * ((j as f64 + 0.5) * k as f64 * pi_over_n).cos())
                    .sum();
                // DCT-II normalisation
                let norm = if k == 0 {
                    (1.0 / n as f64).sqrt()
                } else {
                    (2.0 / n as f64).sqrt()
                };
                sum * norm
            })
            .collect()
    }

    /// Compute the inverse DCT-II (i.e. DCT-III) transform of `coeffs`.
    ///
    /// Reconstructs the original signal from its DCT coefficients.
    pub fn inverse(&self, coeffs: &[f64]) -> Vec<f64> {
        let n = self.n.min(coeffs.len());
        let pi_over_n = std::f64::consts::PI / n as f64;
        let norm_rest = (2.0 / n as f64).sqrt();
        let norm0 = (1.0 / n as f64).sqrt();
        (0..n)
            .map(|j| {
                norm0 * coeffs[0]
                    + coeffs[1..n]
                        .iter()
                        .enumerate()
                        .map(|(ki, &ck)| {
                            let k = ki + 1;
                            norm_rest * ck * ((j as f64 + 0.5) * k as f64 * pi_over_n).cos()
                        })
                        .sum::<f64>()
            })
            .collect()
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
pub struct RandomMeasurementMatrix {
    /// Number of measurements (rows).
    pub m: usize,
    /// Signal length (columns).
    pub n: usize,
    /// Underlying matrix entries stored row-major.
    pub matrix: Vec<Vec<f64>>,
}

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
pub struct BasisPursuit;

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
pub struct OrthogonalMatchingPursuit {
    /// Maximum sparsity (number of non-zero coefficients to recover).
    pub max_k: usize,
}

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
                let corr: f64 = a
                    .iter()
                    .zip(residual.iter())
                    .map(|(row, &ri)| row[j] * ri)
                    .sum::<f64>()
                    .abs();
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
                    ata[si][sj] = a.iter().map(|row| row[ci] * row[cj]).sum();
                }
                atb[si] = a.iter().zip(b.iter()).map(|(row, &bi)| row[ci] * bi).sum();
            }

            // Solve s×s system via Gaussian elimination
            let coeffs = gauss_solve(&ata, &atb);

            // Update x
            for xj in x.iter_mut() {
                *xj = 0.0;
            }
            for (si, &ci) in support.iter().enumerate() {
                x[ci] = coeffs[si];
            }

            // Update residual: r = b - A x
            residual = a
                .iter()
                .zip(b.iter())
                .map(|(row, &bi)| {
                    let ax_i: f64 = row.iter().zip(x.iter()).map(|(aij, xj)| aij * xj).sum();
                    bi - ax_i
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
                let corr: f64 = a
                    .iter()
                    .zip(residual.iter())
                    .map(|(row, &ri)| row[j] * ri)
                    .sum::<f64>()
                    .abs();
                if corr > best_corr {
                    best_corr = corr;
                    best_idx = j;
                }
            }
            support.push(best_idx);
            // Quick residual update (orthogonal projection onto selected atom)
            let col_norm_sq: f64 = a.iter().map(|row| row[best_idx].powi(2)).sum();
            if col_norm_sq < 1e-14 {
                break;
            }
            let proj: f64 = a
                .iter()
                .zip(residual.iter())
                .map(|(row, &ri)| row[best_idx] * ri)
                .sum::<f64>()
                / col_norm_sq;
            for (ri, row) in residual.iter_mut().zip(a.iter()) {
                *ri -= proj * row[best_idx];
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
            let col_slice: Vec<f64> = mat[col][col..n].to_vec();
            for (cell, &cv) in mat[row][col..n].iter_mut().zip(col_slice.iter()) {
                *cell -= factor * cv;
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
pub struct SparsityMetrics;

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
pub struct RecoveryGuarantee;

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

/// Default seed of [`KSvd`] (the ASCII bytes of "KSVD").
const KSVD_DEFAULT_SEED: u64 = 0x4B53_5644;

/// Power iterations per atom in the rank-1 dictionary update.
const KSVD_POWER_ITERATIONS: usize = 10;

/// K-SVD dictionary learning algorithm.
///
/// Alternates between a sparse coding stage (OMP) and a dictionary update stage
/// (rank-1 update of every atom *together with* the coefficients that use it) to
/// learn a possibly overcomplete dictionary adapted to training data.
///
/// Training is deterministic: the initial atoms and the replacements of unused
/// ("dead") atoms are drawn from a generator seeded with [`KSvd::seed`].
///
/// Reference: Aharon, Elad & Bruckstein (2006).
#[derive(Debug, Clone)]
pub struct KSvd {
    /// Number of dictionary atoms (columns).  May exceed the signal dimension
    /// (overcomplete dictionary).
    pub n_atoms: usize,
    /// Target sparsity per signal.
    pub sparsity: usize,
    /// Number of training iterations.
    pub n_iter: usize,
    /// Seed of the pseudo-random generator used to pick the initial atoms and to
    /// re-initialise atoms that no signal uses.
    pub seed: u64,
}

/// Representation errors ‖Y − D C‖_F recorded for one K-SVD iteration by
/// [`KSvd::fit_with_history`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct KSvdIterationError {
    /// Error after the sparse-coding (OMP) stage, with the dictionary of the
    /// previous iteration.
    pub after_coding: f64,
    /// Error after the dictionary-update stage, i.e. after every atom and the
    /// coefficients using it were refreshed.  The update stage never increases
    /// the error, so this is at most `after_coding` (up to rounding).
    pub after_update: f64,
}

impl KSvd {
    /// Create a new K-SVD learner with the default seed.
    ///
    /// - `n_atoms` — number of dictionary atoms
    /// - `sparsity` — maximum number of atoms per signal
    /// - `n_iter` — number of alternating optimisation iterations
    pub fn new(n_atoms: usize, sparsity: usize, n_iter: usize) -> Self {
        Self {
            n_atoms,
            sparsity,
            n_iter,
            seed: KSVD_DEFAULT_SEED,
        }
    }

    /// Return the learner with the random seed replaced by `seed`.
    pub fn with_seed(mut self, seed: u64) -> Self {
        self.seed = seed;
        self
    }

    /// Learn a dictionary from training signals.
    ///
    /// - `signals` — list of training signals (each of the same length `d`)
    ///
    /// Returns the dictionary as `n_atoms` atoms of length `d`: `dict[k]` is the
    /// `k`-th atom (unit ℓ₂ norm) and `dict[k][i]` its `i`-th component.  The
    /// dictionary may be overcomplete (`n_atoms > d`).  Returns an empty vector
    /// if there are no signals or if the signals are empty or differ in length.
    pub fn fit(&self, signals: &[Vec<f64>]) -> Vec<Vec<f64>> {
        self.fit_impl(signals, false).0
    }

    /// Learn a dictionary like [`KSvd::fit`] and also report the representation
    /// error ‖Y − D C‖_F after the coding and after the update stage of every
    /// iteration, which allows monitoring convergence.
    pub fn fit_with_history(
        &self,
        signals: &[Vec<f64>],
    ) -> (Vec<Vec<f64>>, Vec<KSvdIterationError>) {
        self.fit_impl(signals, true)
    }

    /// Shared training loop; `track_errors` enables the (costly) error history.
    fn fit_impl(
        &self,
        signals: &[Vec<f64>],
        track_errors: bool,
    ) -> (Vec<Vec<f64>>, Vec<KSvdIterationError>) {
        let mut history: Vec<KSvdIterationError> = Vec::new();
        let Some(first) = signals.first() else {
            return (Vec::new(), history);
        };
        let d = first.len();
        if d == 0 || signals.iter().any(|y| y.len() != d) {
            return (Vec::new(), history);
        }
        let n_signals = signals.len();
        let n_atoms = self.n_atoms;
        let mut rng = StdRng::seed_from_u64(self.seed);

        // Initialise the dictionary with distinct, randomly chosen training
        // signals (Fisher–Yates shuffle of the signal indices), normalised.
        let mut order: Vec<usize> = (0..n_signals).collect();
        for i in (1..n_signals).rev() {
            let j = rng.random_range(0..=i);
            order.swap(i, j);
        }
        let mut dict: Vec<Vec<f64>> = (0..n_atoms)
            .map(|k| unit_atom(&signals[order[k % n_signals]]))
            .collect();

        let omp = OrthogonalMatchingPursuit::new(self.sparsity);

        for _ in 0..self.n_iter {
            // --- Sparse coding stage: encode each signal with OMP ---
            // We want a sparse code c with D c ≈ y, i.e. A = D (d × n_atoms) and
            // b = y (length d); OMP takes A row-major, so pass D transposed
            // relative to the stored atoms.
            let dict_t: Vec<Vec<f64>> = (0..d)
                .map(|i| dict.iter().map(|atom| atom[i]).collect())
                .collect();
            let mut codes: Vec<Vec<f64>> = signals.iter().map(|y| omp.solve(&dict_t, y)).collect();
            let after_coding = if track_errors {
                representation_error(signals, &dict, &codes)
            } else {
                0.0
            };

            // --- Dictionary update stage: rank-1 update of each atom and of
            //     the coefficients that use it ---
            for k in 0..n_atoms {
                if !ksvd_update_atom(signals, &mut dict, &mut codes, k) {
                    // Re-initialise a dead atom with a random training signal.
                    let pick = rng.random_range(0..n_signals);
                    dict[k] = unit_atom(&signals[pick]);
                }
            }

            if track_errors {
                history.push(KSvdIterationError {
                    after_coding,
                    after_update: representation_error(signals, &dict, &codes),
                });
            }
        }

        (dict, history)
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

/// Copy of `signal` scaled to unit ℓ₂ norm (a zero signal stays zero).
fn unit_atom(signal: &[f64]) -> Vec<f64> {
    let norm = l2_norm(signal);
    if norm > 0.0 && norm.is_finite() {
        signal.iter().map(|v| v / norm).collect()
    } else {
        signal.to_vec()
    }
}

/// Rank-1 K-SVD update of atom `k` and of the coefficients that use it.
///
/// With `E_k = Y_k − Σ_{j≠k} d_j c_jᵀ` the residual of the signals `Y_k` that use
/// atom `k` (computed from the *current* atoms and codes), the best rank-1 fit
/// `E_k ≈ d_k c_kᵀ` is obtained by taking `d_k` as the dominant direction of
/// `E_k` and refreshing the coefficients to the projections `c_k = E_k d_k`.
/// Refreshing the coefficients together with the atom is essential: the other
/// atoms' residuals are built from `codes`, and stale coefficients would also
/// make the stage non-monotone.
///
/// Returns `false` (leaving everything untouched) if no signal uses atom `k`.
fn ksvd_update_atom(
    signals: &[Vec<f64>],
    dict: &mut [Vec<f64>],
    codes: &mut [Vec<f64>],
    k: usize,
) -> bool {
    // Signals that use atom k.
    let using: Vec<usize> = (0..signals.len())
        .filter(|&s| codes[s][k].abs() > 1e-14)
        .collect();
    if using.is_empty() {
        return false;
    }

    let e_rows: Vec<Vec<f64>> = using
        .iter()
        .map(|&s| {
            let mut e = signals[s].clone();
            for (j, dict_j) in dict.iter().enumerate() {
                let coef = codes[s][j];
                if j == k || coef == 0.0 {
                    continue;
                }
                for (ei, &dji) in e.iter_mut().zip(dict_j.iter()) {
                    *ei -= coef * dji;
                }
            }
            e
        })
        .collect();

    // The power iteration is warm-started at the current (unit-norm) atom, so
    // the rank-1 fit never gets worse than with the previous atom.
    let atom = dominant_direction(&e_rows, &dict[k]);
    for (row, &s) in e_rows.iter().zip(using.iter()) {
        codes[s][k] = row.iter().zip(atom.iter()).map(|(a, b)| a * b).sum();
    }
    dict[k] = atom;
    true
}

/// Dominant right singular direction of the matrix whose rows are `rows`.
///
/// Runs power iteration on `RᵀR` starting from the unit vector `start`.  The
/// Rayleigh quotient `Σ_r (row_r · v)²` never decreases along the iteration, so
/// the result fits `rows` at least as well (as a rank-1 direction) as `start`.
/// If every row is orthogonal to the current iterate the iteration stops and
/// keeps it.
fn dominant_direction(rows: &[Vec<f64>], start: &[f64]) -> Vec<f64> {
    let mut atom = start.to_vec();
    for _ in 0..KSVD_POWER_ITERATIONS {
        let mut next = vec![0.0_f64; atom.len()];
        for row in rows {
            let projection: f64 = row.iter().zip(atom.iter()).map(|(a, b)| a * b).sum();
            for (n_i, &r_i) in next.iter_mut().zip(row.iter()) {
                *n_i += r_i * projection;
            }
        }
        let norm = l2_norm(&next);
        if norm == 0.0 || !norm.is_finite() {
            break;
        }
        for v in next.iter_mut() {
            *v /= norm;
        }
        atom = next;
    }
    atom
}

/// Frobenius norm ‖Y − D C‖_F of the representation residual.
///
/// `dict[k]` is atom `k` and `codes[s][k]` the coefficient of atom `k` in the
/// code of signal `s`.
fn representation_error(signals: &[Vec<f64>], dict: &[Vec<f64>], codes: &[Vec<f64>]) -> f64 {
    let mut total = 0.0_f64;
    for (signal, code) in signals.iter().zip(codes.iter()) {
        let mut residual = signal.clone();
        for (atom, &coef) in dict.iter().zip(code.iter()) {
            if coef == 0.0 {
                continue;
            }
            for (r, &a) in residual.iter_mut().zip(atom.iter()) {
                *r -= coef * a;
            }
        }
        total += residual.iter().map(|r| r * r).sum::<f64>();
    }
    total.sqrt()
}

// ─────────────────────────────────────────────────────────────────────────────
// MRI-like Compressed Sensing
// ─────────────────────────────────────────────────────────────────────────────

/// Compressed sensing reconstruction for MRI-like k-space data.
///
/// In MRI, measurements are taken in Fourier (k-space) domain.  This module
/// provides a simplified 1-D model: the signal is sparse in the DCT domain,
/// and k-space samples are random Fourier measurements.
pub struct MriCompressedSensing {
    /// Signal length.
    pub n: usize,
    /// Number of k-space samples (measurements).
    pub m: usize,
}

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
pub struct SparseSignal;

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
#[path = "compressed_sensing_tests.rs"]
mod tests;
