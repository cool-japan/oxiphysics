// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Eigenvalue analysis for FEM structural dynamics.
//!
//! Provides multiple eigensolver strategies for the generalized eigenvalue
//! problem K·φ = λM·φ arising in structural dynamics, plus utilities for
//! modal analysis, frequency response, and model-order reduction.
//!
//! # Included solvers
//! - [`PowerIteration`] — dominant eigenvalue via power method with deflation
//! - [`InverseIteration`] — eigenvalue nearest a shift via inverse iteration
//!   (pivoted LU)
//! - [`LanczosMethod`] — Lanczos with full reorthogonalization for the k
//!   smallest eigenvalues of a symmetric matrix
//! - [`GeneralizedEigenvalue`] — K·φ = λM·φ via Cholesky of M and a dense
//!   symmetric eigendecomposition
//! - [`ModalAnalysis`] — natural frequencies, mode shapes, modal participation
//! - [`FrequencyResponse`] — steady-state harmonic response with hysteretic
//!   damping
//! - [`RayleighQuotient`] — upper bound on fundamental frequency
//! - [`EigensolverConfig`] — iteration / convergence configuration
//! - [`SpectralNorm`] — largest singular value via power iteration
//! - [`ModalTruncation`] — reduced-order model from retained modes
//! - [`lobpcg_solve`] — LOBPCG with optional shift-invert for sparse `CsrMatrix` pencils

pub mod lobpcg;
pub use lobpcg::{EigensolverError, LobpcgConfig, LobpcgResult, lobpcg_solve};

// ─────────────────────────────────────────────────────────────────────────────
// Utilities
// ─────────────────────────────────────────────────────────────────────────────

/// Dot product of two slices.
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

/// 2-norm of a slice.
fn norm2(v: &[f64]) -> f64 {
    dot(v, v).sqrt()
}

/// Normalize a vector in-place; returns the previous norm.
fn normalize(v: &mut [f64]) -> f64 {
    let n = norm2(v);
    if n > 1e-300 {
        for x in v.iter_mut() {
            *x /= n;
        }
    }
    n
}

/// Dense matrix-vector product y = A x  where A is stored row-major (n×n).
fn dense_matvec(a: &[f64], x: &[f64], n: usize) -> Vec<f64> {
    let mut y = vec![0.0_f64; n];
    for i in 0..n {
        for j in 0..n {
            y[i] += a[i * n + j] * x[j];
        }
    }
    y
}

/// Largest absolute row sum (the infinity norm) of a dense row-major `n × n`
/// matrix.  For a symmetric matrix this bounds the spectral radius from above.
fn infinity_norm(a: &[f64], n: usize) -> f64 {
    a.chunks_exact(n.max(1))
        .map(|row| row.iter().map(|v| v.abs()).sum::<f64>())
        .fold(0.0_f64, f64::max)
}

/// Maximum number of QL iterations spent on one eigenvalue of a tridiagonal
/// matrix before the routine gives up.
const TRIDIAGONAL_MAX_ITER: usize = 60;

/// Eigenvalues of a symmetric tridiagonal matrix together with the last
/// component of each eigenvector.
///
/// `diag` is the diagonal (length `n`) and `off` the sub-diagonal (length
/// `n − 1`).  The implicit QL algorithm with Wilkinson-type shifts (EISPACK
/// `tql2`) is run while only the last row of the accumulated rotations is
/// tracked, which costs `O(n²)` instead of the `O(n³)` of a full eigenvector
/// computation.  The last components are exactly what the residual of a Ritz
/// pair needs: for a Lanczos step with residual coupling `β`, the residual norm
/// of the Ritz pair `(θᵢ, sᵢ)` is `β·|sᵢ[n−1]|`.
///
/// Returns `(eigenvalues ascending, matching last eigenvector components)`, or
/// `None` if an eigenvalue does not converge within
/// [`TRIDIAGONAL_MAX_ITER`] iterations.
fn tridiagonal_eigen(diag: &[f64], off: &[f64]) -> Option<(Vec<f64>, Vec<f64>)> {
    let n = diag.len();
    if n == 0 {
        return Some((Vec::new(), Vec::new()));
    }
    let mut d = diag.to_vec();
    // e[i] couples rows i and i + 1; e[n - 1] is a zero sentinel.
    let mut e = vec![0.0_f64; n];
    for (slot, &value) in e.iter_mut().zip(off.iter().take(n - 1)) {
        *slot = value;
    }
    // Last row of the eigenvector matrix; starts as the last row of the identity.
    let mut z = vec![0.0_f64; n];
    z[n - 1] = 1.0;

    let mut shift_total = 0.0_f64;
    let mut tst1 = 0.0_f64;
    for l in 0..n {
        tst1 = tst1.max(d[l].abs() + e[l].abs());
        let mut m = l;
        while m < n - 1 && e[m].abs() > f64::EPSILON * tst1 {
            m += 1;
        }
        if m > l {
            let mut iterations = 0usize;
            loop {
                iterations += 1;
                if iterations > TRIDIAGONAL_MAX_ITER {
                    return None;
                }
                // Implicit shift.
                let mut g = d[l];
                let mut p = (d[l + 1] - g) / (2.0 * e[l]);
                let mut r = p.hypot(1.0);
                if p < 0.0 {
                    r = -r;
                }
                d[l] = e[l] / (p + r);
                d[l + 1] = e[l] * (p + r);
                let dl1 = d[l + 1];
                let mut h = g - d[l];
                for value in d.iter_mut().skip(l + 2) {
                    *value -= h;
                }
                shift_total += h;

                // Implicit QL transformation.
                p = d[m];
                let mut c = 1.0_f64;
                let mut c2 = c;
                let mut c3 = c;
                let el1 = e[l + 1];
                let mut s = 0.0_f64;
                let mut s2 = 0.0_f64;
                for i in (l..m).rev() {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c * e[i];
                    h = c * p;
                    r = p.hypot(e[i]);
                    e[i + 1] = s * r;
                    s = e[i] / r;
                    c = p / r;
                    p = c * d[i] - s * g;
                    d[i + 1] = h + s * (c * g + s * d[i]);
                    let z_next = z[i + 1];
                    z[i + 1] = s * z[i] + c * z_next;
                    z[i] = c * z[i] - s * z_next;
                }
                p = -s * s2 * c3 * el1 * e[l] / dl1;
                e[l] = s * p;
                d[l] = c * p;
                if e[l].abs() <= f64::EPSILON * tst1 {
                    break;
                }
            }
        }
        d[l] += shift_total;
        e[l] = 0.0;
    }
    if d.iter().chain(z.iter()).any(|v| !v.is_finite()) {
        return None;
    }
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| d[a].total_cmp(&d[b]));
    let values = order.iter().map(|&i| d[i]).collect();
    let last = order.iter().map(|&i| z[i]).collect();
    Some((values, last))
}

/// Ritz values (ascending) and the last components of the matching Ritz vectors
/// of the symmetric tridiagonal matrix with diagonal `alpha` and sub-diagonal
/// `beta`, computed with the dense Jacobi method.  Slower than
/// [`tridiagonal_eigen`] but unconditionally convergent; used as its fallback.
fn ritz_pairs_jacobi(alpha: &[f64], beta: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let m = alpha.len();
    if m == 0 {
        return (Vec::new(), Vec::new());
    }
    let scale = alpha
        .iter()
        .chain(beta.iter())
        .fold(0.0_f64, |s, v| s.max(v.abs()));
    let scale = if scale > 0.0 { scale } else { 1.0 };
    let mut t = vec![vec![0.0_f64; m]; m];
    for (i, &a) in alpha.iter().enumerate() {
        t[i][i] = a / scale;
    }
    for (i, &b) in beta.iter().enumerate().take(m - 1) {
        t[i][i + 1] = b / scale;
        t[i + 1][i] = b / scale;
    }
    let (values, vectors) = crate::eigen::jacobi_eigen_dense(&t, m);
    let last = vectors
        .iter()
        .map(|v| v.last().copied().unwrap_or(0.0))
        .collect();
    (values.iter().map(|v| v * scale).collect(), last)
}

/// Ritz values (ascending) and last Ritz-vector components of a Lanczos
/// tridiagonal matrix (see [`tridiagonal_eigen`]).
fn ritz_pairs(alpha: &[f64], beta: &[f64]) -> (Vec<f64>, Vec<f64>) {
    tridiagonal_eigen(alpha, beta).unwrap_or_else(|| ritz_pairs_jacobi(alpha, beta))
}

/// Deterministic pseudo-random vector with entries in `[-1, 1)` (xorshift64*),
/// used for the Lanczos start and restart vectors.
fn pseudo_random_vector(n: usize, seed: u64) -> Vec<f64> {
    let mut state = seed | 1;
    (0..n)
        .map(|_| {
            state ^= state >> 12;
            state ^= state << 25;
            state ^= state >> 27;
            let bits = state.wrapping_mul(0x2545_F491_4F6C_DD1D) >> 11;
            (bits as f64) / ((1u64 << 52) as f64) - 1.0
        })
        .collect()
}

/// Remove from `w` its components along every (orthonormal) vector of `basis`.
/// Gram–Schmidt is applied twice so that the result is orthogonal to working
/// precision even when `w` is nearly in the span of the basis.
fn orthogonalize_against(w: &mut [f64], basis: &[Vec<f64>]) {
    for _ in 0..2 {
        for q in basis {
            let h = dot(w, q);
            for (wi, qi) in w.iter_mut().zip(q) {
                *wi -= h * qi;
            }
        }
    }
}

/// A fresh unit vector orthogonal to `basis`, or `None` if `basis` already
/// spans the whole `n`-dimensional space.
fn restart_vector(n: usize, basis: &[Vec<f64>], seed: &mut u64) -> Option<Vec<f64>> {
    if basis.len() >= n {
        return None;
    }
    for _ in 0..8 {
        *seed = seed.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut v = pseudo_random_vector(n, *seed);
        let before = norm2(&v);
        orthogonalize_against(&mut v, basis);
        let after = norm2(&v);
        if after > 1e-8 * before {
            for x in v.iter_mut() {
                *x /= after;
            }
            return Some(v);
        }
    }
    None
}

/// Relative pivot threshold (times `max|A|`) below which the dynamic stiffness
/// matrix of the frequency-response solver is declared singular.
const FRF_PIVOT_TOL: f64 = 1e-14;

/// In-place dense LU factorization with partial (row) pivoting of an `n × n`
/// row-major matrix: `P·A = L·U`.
///
/// `L` (unit lower triangle) and `U` overwrite `a`.  The returned vector is the
/// row permutation: row `i` of `P·A` is row `perm[i]` of `A`.  `None` means the
/// matrix is singular to working precision: it holds a non-finite entry, or a
/// pivot is not larger than `rel_tol · max|A|` (with `rel_tol = 0` only an exact
/// zero pivot is rejected).
fn lu_factor(a: &mut [f64], n: usize, rel_tol: f64) -> Option<Vec<usize>> {
    if a.iter().any(|v| !v.is_finite()) {
        return None;
    }
    let threshold = rel_tol * a.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    let mut perm: Vec<usize> = (0..n).collect();
    for k in 0..n {
        let mut pivot_row = k;
        let mut pivot_abs = a[k * n + k].abs();
        for i in (k + 1)..n {
            let candidate = a[i * n + k].abs();
            if candidate > pivot_abs {
                pivot_abs = candidate;
                pivot_row = i;
            }
        }
        if pivot_abs <= threshold {
            return None;
        }
        if pivot_row != k {
            for j in 0..n {
                a.swap(k * n + j, pivot_row * n + j);
            }
            perm.swap(k, pivot_row);
        }
        let pivot = a[k * n + k];
        for i in (k + 1)..n {
            let factor = a[i * n + k] / pivot;
            a[i * n + k] = factor;
            for j in (k + 1)..n {
                let upper = a[k * n + j];
                a[i * n + j] -= factor * upper;
            }
        }
    }
    Some(perm)
}

/// Solve `A·x = b` from the factors and permutation returned by [`lu_factor`].
fn lu_solve(lu: &[f64], perm: &[usize], b: &[f64], n: usize) -> Vec<f64> {
    let mut x: Vec<f64> = perm.iter().map(|&p| b[p]).collect();
    // Forward substitution L·y = P·b (unit lower triangle).
    for i in 1..n {
        for j in 0..i {
            x[i] -= lu[i * n + j] * x[j];
        }
    }
    // Back substitution U·x = y.
    for i in (0..n).rev() {
        for j in (i + 1)..n {
            x[i] -= lu[i * n + j] * x[j];
        }
        x[i] /= lu[i * n + i];
    }
    x
}

// ─────────────────────────────────────────────────────────────────────────────
// EigensolverConfig
// ─────────────────────────────────────────────────────────────────────────────

/// Configuration shared by all eigensolvers.
#[derive(Debug, Clone)]
pub struct EigensolverConfig {
    /// Maximum number of outer iterations.
    pub max_iter: usize,
    /// Convergence tolerance on the residual or eigenvalue change.
    pub tolerance: f64,
    /// Spectral shift σ used in shift-and-invert strategies.
    pub shift: f64,
    /// Number of eigenvalues to compute.
    pub n_modes: usize,
}

impl Default for EigensolverConfig {
    fn default() -> Self {
        Self {
            max_iter: 500,
            tolerance: 1e-10,
            shift: 0.0,
            n_modes: 5,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PowerIteration
// ─────────────────────────────────────────────────────────────────────────────

/// Power iteration for the dominant eigenvalue/eigenvector of a dense matrix.
///
/// Deflation is used to extract subsequent eigenvalues after the dominant one.
/// Storage: the matrix `a` is row-major, size n×n (flat `Vec`f64`).
pub struct PowerIteration {
    /// Row-major n×n matrix.
    pub a: Vec<f64>,
    /// Matrix dimension.
    pub n: usize,
    /// Solver configuration.
    pub config: EigensolverConfig,
}

impl PowerIteration {
    /// Create a new power iteration solver.
    pub fn new(a: Vec<f64>, n: usize, config: EigensolverConfig) -> Self {
        assert_eq!(a.len(), n * n, "matrix must be n×n");
        Self { a, n, config }
    }

    /// Compute the dominant (largest-magnitude) eigenvalue and eigenvector.
    ///
    /// Returns `(eigenvalue, eigenvector)`.
    pub fn dominant(&self) -> (f64, Vec<f64>) {
        let n = self.n;
        // Use a non-axis-aligned starting vector so it has components in all eigenspaces.
        let scale = (n as f64).sqrt().recip();
        let mut v: Vec<f64> = (0..n).map(|i| (i + 1) as f64 * scale).collect();
        normalize(&mut v);
        let mut lambda = 0.0_f64;
        for _ in 0..self.config.max_iter {
            let w = dense_matvec(&self.a, &v, n);
            let new_lambda = dot(&v, &w);
            let mut new_v = w;
            let nrm = normalize(&mut new_v);
            if nrm < 1e-300 {
                break;
            }
            if (new_lambda - lambda).abs() < self.config.tolerance {
                return (new_lambda, new_v);
            }
            lambda = new_lambda;
            v = new_v;
        }
        (lambda, v)
    }

    /// Compute the `k` largest eigenvalues via power iteration with deflation.
    ///
    /// Returns vectors of `(eigenvalue, eigenvector)` pairs sorted descending.
    pub fn deflated_eigenvalues(&self, k: usize) -> Vec<(f64, Vec<f64>)> {
        let n = self.n;
        let mut a_deflated = self.a.clone();
        let mut results = Vec::with_capacity(k);
        for _ in 0..k {
            let solver = PowerIteration {
                a: a_deflated.clone(),
                n,
                config: self.config.clone(),
            };
            let (lam, vec) = solver.dominant();
            // Deflate: A ← A - λ v vᵀ
            for i in 0..n {
                for j in 0..n {
                    a_deflated[i * n + j] -= lam * vec[i] * vec[j];
                }
            }
            results.push((lam, vec));
        }
        results
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// InverseIteration
// ─────────────────────────────────────────────────────────────────────────────

/// Inverse iteration for the smallest (or nearest-to-shift) eigenvalue.
///
/// Uses LU factorization with partial pivoting of (A − σ·I) and iterates with
/// back-substitution.
pub struct InverseIteration {
    /// Row-major n×n stiffness matrix.
    pub a: Vec<f64>,
    /// Matrix dimension.
    pub n: usize,
    /// Solver configuration (shift σ is used here).
    pub config: EigensolverConfig,
}

impl InverseIteration {
    /// Create a new inverse-iteration solver.
    pub fn new(a: Vec<f64>, n: usize, config: EigensolverConfig) -> Self {
        assert_eq!(a.len(), n * n);
        Self { a, n, config }
    }

    /// Solve for the eigenvalue nearest to `config.shift`.
    ///
    /// Returns `(eigenvalue, eigenvector)`.  Uses Rayleigh quotient to
    /// estimate the converged eigenvalue.
    pub fn solve(&self) -> Result<(f64, Vec<f64>), String> {
        let n = self.n;
        let sigma = self.config.shift;
        // Form shifted matrix (A - σI)
        let mut lu = self.a.clone();
        for i in 0..n {
            lu[i * n + i] -= sigma;
        }
        // Only an exactly singular shifted matrix is rejected: a shift that is
        // merely very close to an eigenvalue is what makes the iteration fast.
        let perm =
            lu_factor(&mut lu, n, 0.0).ok_or_else(|| "singular shifted matrix".to_string())?;

        // Non-axis-aligned starting vector
        let scale = (n as f64).sqrt().recip();
        let mut v: Vec<f64> = (0..n).map(|i| (i + 1) as f64 * scale).collect();
        normalize(&mut v);

        for _ in 0..self.config.max_iter {
            let w = lu_solve(&lu, &perm, &v, n);
            let mut new_v = w;
            normalize(&mut new_v);
            // Check convergence: angle between successive iterates
            let angle_cos = dot(&new_v, &v).abs();
            let converged = (1.0 - angle_cos) < self.config.tolerance;
            v = new_v;
            if converged {
                break;
            }
        }
        // Rayleigh quotient: λ = vᵀ A v / vᵀ v (vᵀv = 1 since normalized)
        let av = dense_matvec(&self.a, &v, n);
        let eigenvalue = dot(&v, &av);
        Ok((eigenvalue, v))
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// LanczosMethod
// ─────────────────────────────────────────────────────────────────────────────

/// Breakdown threshold of the Lanczos recurrence, relative to `‖A‖∞`: a
/// residual vector shorter than this means the Krylov space is invariant.
const LANCZOS_BREAKDOWN_TOL: f64 = 1e-13;

/// Floor of the Ritz residual test, relative to `‖A‖∞`: residuals cannot be
/// pushed below the rounding level of the matrix-vector products.
const LANCZOS_RESIDUAL_FLOOR: f64 = 1e-14;

/// Seed of the deterministic pseudo-random Lanczos start vector.
const LANCZOS_SEED: u64 = 0x9E37_79B9_7F4A_7C15;

/// Lanczos tridiagonalization for extracting the k smallest eigenvalues of a
/// symmetric matrix.
///
/// Builds a Krylov basis with **full reorthogonalization** and extracts the Ritz
/// values from the resulting tridiagonal matrix.  The basis is grown from
/// `min(n, max(2k + 20, 3k))` vectors until the `k` smallest Ritz values have
/// converged (Ritz residual `β·|s_last| ≤ tolerance·|θ|`) or the basis spans the
/// whole space, in which case the Ritz values are the exact eigenvalues.
///
/// A single-vector Krylov method resolves each *distinct* eigenvalue once, so a
/// multiple eigenvalue is reported with its full multiplicity only when the
/// Krylov space closes early and the recurrence restarts, or when the basis
/// reaches dimension `n`.  Use [`GeneralizedEigenvalue`] when exact
/// multiplicities of the lowest modes are required.
pub struct LanczosMethod {
    /// Row-major n×n symmetric matrix.
    pub a: Vec<f64>,
    /// Matrix dimension.
    pub n: usize,
    /// Solver configuration (`n_modes` determines k, `tolerance` the Ritz
    /// residual tolerance; the basis size is bounded by `n`, not by `max_iter`).
    pub config: EigensolverConfig,
}

impl LanczosMethod {
    /// Create a new Lanczos solver.
    pub fn new(a: Vec<f64>, n: usize, config: EigensolverConfig) -> Self {
        assert_eq!(a.len(), n * n);
        Self { a, n, config }
    }

    /// Run Lanczos iteration and return `min(config.n_modes, n)` eigenvalues.
    ///
    /// Returns the smallest eigenvalues sorted ascending.  A matrix with
    /// non-finite entries yields NaN values.
    pub fn eigenvalues(&self) -> Vec<f64> {
        let n = self.n;
        let k = self.config.n_modes.min(n);
        if k == 0 {
            return Vec::new();
        }
        if self.a.iter().any(|v| !v.is_finite()) {
            return vec![f64::NAN; k];
        }
        let a_norm = infinity_norm(&self.a, n);
        if a_norm == 0.0 {
            return vec![0.0; k];
        }
        let breakdown_tol = LANCZOS_BREAKDOWN_TOL * a_norm;
        let mut seed = LANCZOS_SEED;
        let mut basis: Vec<Vec<f64>> = Vec::new();
        let Some(start) = restart_vector(n, &basis, &mut seed) else {
            return vec![f64::NAN; k];
        };
        basis.push(start);
        let mut alpha: Vec<f64> = Vec::new();
        // beta[j] couples q_j and q_{j+1}; it is zero where the Krylov space
        // closed and the recurrence was restarted with a fresh vector.
        let mut beta: Vec<f64> = Vec::new();
        let mut m_limit = n.min((2 * k + 20).max(3 * k));
        let mut exhausted = false;

        loop {
            while alpha.len() < m_limit {
                let j = alpha.len();
                let Some(q_j) = basis.get(j) else {
                    exhausted = true;
                    break;
                };
                let mut w = dense_matvec(&self.a, q_j, n);
                let a_j = dot(&w, q_j);
                alpha.push(a_j);
                for (wi, qi) in w.iter_mut().zip(q_j) {
                    *wi -= a_j * qi;
                }
                if j > 0 && beta[j - 1] != 0.0 {
                    let b_prev = beta[j - 1];
                    for (wi, qi) in w.iter_mut().zip(&basis[j - 1]) {
                        *wi -= b_prev * qi;
                    }
                }
                orthogonalize_against(&mut w, &basis);
                let b_j = norm2(&w);
                if b_j > breakdown_tol {
                    for wi in w.iter_mut() {
                        *wi /= b_j;
                    }
                    beta.push(b_j);
                    basis.push(w);
                } else {
                    // Invariant subspace: close the block and restart with a
                    // vector orthogonal to everything found so far.
                    beta.push(0.0);
                    if let Some(next) = restart_vector(n, &basis, &mut seed) {
                        basis.push(next);
                    }
                }
            }

            let m = alpha.len();
            let (ritz, last) = ritz_pairs(&alpha, &beta[..m.saturating_sub(1)]);
            let coupling = beta.get(m.saturating_sub(1)).copied().unwrap_or(0.0);
            let converged = coupling > 0.0
                && ritz.iter().zip(&last).take(k).all(|(&theta, &s_last)| {
                    let tolerance =
                        (self.config.tolerance * theta.abs()).max(LANCZOS_RESIDUAL_FLOOR * a_norm);
                    coupling * s_last.abs() <= tolerance
                });
            if m >= n || exhausted || converged {
                return ritz.into_iter().take(k).collect();
            }
            m_limit = n.min(m + (m / 2).max(10));
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// GeneralizedEigenvalue
// ─────────────────────────────────────────────────────────────────────────────

/// Relative asymmetry (times `max|A|`) tolerated in a stiffness or mass matrix
/// before [`GeneralizedEigenvalue::solve_smallest`] rejects it.
const SYMMETRY_TOL: f64 = 1e-8;

/// Generalized eigenvalue solver for K·φ = λM·φ.
///
/// Uses Cholesky factorization of M to convert to a standard eigenvalue problem:
/// `L⁻¹ K L⁻ᵀ ỹ = λ ỹ`, where M = L Lᵀ, and solves that problem with a dense
/// symmetric eigendecomposition (`O(n³)`), so every requested mode is exact and
/// multiple eigenvalues keep their multiplicity.
pub struct GeneralizedEigenvalue {
    /// Row-major n×n stiffness matrix K.
    pub k_mat: Vec<f64>,
    /// Row-major n×n positive-definite mass matrix M.
    pub m_mat: Vec<f64>,
    /// Matrix dimension.
    pub n: usize,
    /// Solver configuration.
    pub config: EigensolverConfig,
}

impl GeneralizedEigenvalue {
    /// Create a new generalized eigenvalue solver.
    pub fn new(k_mat: Vec<f64>, m_mat: Vec<f64>, n: usize, config: EigensolverConfig) -> Self {
        assert_eq!(k_mat.len(), n * n);
        assert_eq!(m_mat.len(), n * n);
        Self {
            k_mat,
            m_mat,
            n,
            config,
        }
    }

    /// Cholesky factor L (lower triangular) of a positive-definite matrix.
    fn cholesky(m: &[f64], n: usize) -> Result<Vec<f64>, String> {
        let mut l = vec![0.0_f64; n * n];
        for i in 0..n {
            for j in 0..=i {
                let mut s: f64 = m[i * n + j];
                for p in 0..j {
                    s -= l[i * n + p] * l[j * n + p];
                }
                if i == j {
                    if s <= 0.0 {
                        return Err("mass matrix is not positive definite".into());
                    }
                    l[i * n + j] = s.sqrt();
                } else {
                    l[i * n + j] = s / l[j * n + j];
                }
            }
        }
        Ok(l)
    }

    /// Solve L x = b (forward substitution, L lower triangular).
    fn fwd_sub(l: &[f64], b: &[f64], n: usize) -> Vec<f64> {
        let mut x = b.to_vec();
        for i in 0..n {
            for j in 0..i {
                x[i] -= l[i * n + j] * x[j];
            }
            x[i] /= l[i * n + i];
        }
        x
    }

    /// Solve Lᵀ x = b (back substitution).
    fn bwd_sub(l: &[f64], b: &[f64], n: usize) -> Vec<f64> {
        let mut x = b.to_vec();
        for i in (0..n).rev() {
            for j in (i + 1)..n {
                x[i] -= l[j * n + i] * x[j];
            }
            x[i] /= l[i * n + i];
        }
        x
    }

    /// Check that `a` (row-major `n × n`) is finite and symmetric to within
    /// [`SYMMETRY_TOL`] of its largest entry.
    fn check_symmetric(a: &[f64], n: usize, name: &str) -> Result<(), String> {
        if a.iter().any(|v| !v.is_finite()) {
            return Err(format!("{name} matrix contains non-finite entries"));
        }
        let tolerance = SYMMETRY_TOL * a.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
        for i in 0..n {
            for j in (i + 1)..n {
                if (a[i * n + j] - a[j * n + i]).abs() > tolerance {
                    return Err(format!(
                        "{name} matrix is not symmetric: entries ({i},{j}) and ({j},{i}) differ"
                    ));
                }
            }
        }
        Ok(())
    }

    /// Compute the `config.n_modes` smallest generalized eigenvalues.
    ///
    /// Returns `(eigenvalues, eigenvectors)`: the eigenvalues λ = ω² ascending
    /// and the matching eigenvectors φ of length `n`, normalised to unit modal
    /// mass (`φᵀ M φ = 1`, hence `φᵀ K φ = λ`) with the component of largest
    /// magnitude positive.  The standard problem `C ỹ = λ ỹ`, `C = L⁻¹ K L⁻ᵀ`,
    /// is solved completely by a dense symmetric eigendecomposition and the
    /// eigenvectors are mapped back with `φ = L⁻ᵀ ỹ`; `config.tolerance` and
    /// `config.max_iter` are therefore not used.
    ///
    /// # Errors
    ///
    /// Returns an error if K or M is not finite or not symmetric, or if M is not
    /// positive definite.
    pub fn solve_smallest(&self) -> Result<(Vec<f64>, Vec<Vec<f64>>), String> {
        let n = self.n;
        Self::check_symmetric(&self.k_mat, n, "stiffness")?;
        Self::check_symmetric(&self.m_mat, n, "mass")?;
        let l = Self::cholesky(&self.m_mat, n)?;
        let k = self.config.n_modes.min(n);
        if k == 0 {
            return Ok((Vec::new(), Vec::new()));
        }

        // Form C = L⁻¹ K L⁻ᵀ  (standard symmetric eigenvalue problem)
        // First compute B = K L⁻ᵀ  via column operations
        let mut c = vec![0.0_f64; n * n];
        for j in 0..n {
            // Column j of L⁻ᵀ is e_j solved from Lᵀ x = e_j
            let mut ej = vec![0.0_f64; n];
            ej[j] = 1.0;
            let linvt_col = Self::bwd_sub(&l, &ej, n);
            // B[:,j] = K * linvt_col
            let b_col = dense_matvec(&self.k_mat, &linvt_col, n);
            // C[:,j] = L⁻¹ * b_col
            let c_col = Self::fwd_sub(&l, &b_col, n);
            for i in 0..n {
                c[i * n + j] = c_col[i];
            }
        }
        // C is symmetric in exact arithmetic; remove the rounding asymmetry.
        for i in 0..n {
            for j in (i + 1)..n {
                let mean = 0.5 * (c[i * n + j] + c[j * n + i]);
                c[i * n + j] = mean;
                c[j * n + i] = mean;
            }
        }

        // The Jacobi kernel works with absolute thresholds, so hand it a matrix
        // whose largest entry is one and rescale the eigenvalues afterwards.
        let scale = c.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
        let scale = if scale > 0.0 { scale } else { 1.0 };
        let scaled: Vec<Vec<f64>> = c
            .chunks_exact(n)
            .map(|row| row.iter().map(|v| v / scale).collect())
            .collect();
        let (values, vectors) = crate::eigen::jacobi_eigen_dense(&scaled, n);

        // Back-transform the eigenvectors: φ = L⁻ᵀ ỹ.
        let mut eigenvalues = Vec::with_capacity(k);
        let mut eigenvectors = Vec::with_capacity(k);
        for (value, y) in values.iter().zip(vectors.iter()).take(k) {
            eigenvalues.push(value * scale);
            let mut phi = Self::bwd_sub(&l, y, n);
            let dominant = phi.iter().fold(
                0.0_f64,
                |best, &v| if v.abs() > best.abs() { v } else { best },
            );
            if dominant < 0.0 {
                phi.iter_mut().for_each(|v| *v = -*v);
            }
            eigenvectors.push(phi);
        }
        Ok((eigenvalues, eigenvectors))
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// RayleighQuotient
// ─────────────────────────────────────────────────────────────────────────────

/// Rayleigh quotient: upper bound on the fundamental natural frequency squared.
///
/// ρ(x) = xᵀ K x / xᵀ M x ≥ ω₁²
pub struct RayleighQuotient {
    /// Row-major n×n stiffness matrix K.
    pub k_mat: Vec<f64>,
    /// Row-major n×n mass matrix M.
    pub m_mat: Vec<f64>,
    /// Matrix dimension.
    pub n: usize,
}

impl RayleighQuotient {
    /// Create a Rayleigh quotient evaluator.
    pub fn new(k_mat: Vec<f64>, m_mat: Vec<f64>, n: usize) -> Self {
        Self { k_mat, m_mat, n }
    }

    /// Evaluate ρ(x) = xᵀ K x / xᵀ M x.
    pub fn evaluate(&self, x: &[f64]) -> f64 {
        let kx = dense_matvec(&self.k_mat, x, self.n);
        let mx = dense_matvec(&self.m_mat, x, self.n);
        let num: f64 = x.iter().zip(kx.iter()).map(|(a, b)| a * b).sum();
        let den: f64 = x.iter().zip(mx.iter()).map(|(a, b)| a * b).sum();
        if den.abs() < 1e-300 {
            f64::INFINITY
        } else {
            num / den
        }
    }

    /// Iteratively minimise the Rayleigh quotient using steepest descent.
    ///
    /// Returns the minimised quotient (approximation to ω₁²) and the
    /// corresponding trial vector.
    pub fn minimise(&self, n_iter: usize) -> (f64, Vec<f64>) {
        let n = self.n;
        let mut x: Vec<f64> = (0..n).map(|i| if i == 0 { 1.0 } else { 0.0 }).collect();
        normalize(&mut x);
        let mut rho = self.evaluate(&x);

        for _ in 0..n_iter {
            let kx = dense_matvec(&self.k_mat, &x, n);
            let mx = dense_matvec(&self.m_mat, &x, n);
            // Gradient of ρ: g = 2(Kx - ρ Mx) / xᵀMx
            let xmx: f64 = x.iter().zip(mx.iter()).map(|(a, b)| a * b).sum();
            let grad: Vec<f64> = kx
                .iter()
                .zip(mx.iter())
                .map(|(&ki, &mi)| 2.0 * (ki - rho * mi) / (xmx + 1e-300))
                .collect();
            // Line-search step
            let step = 0.01;
            let mut xnew: Vec<f64> = x
                .iter()
                .zip(grad.iter())
                .map(|(xi, gi)| xi - step * gi)
                .collect();
            normalize(&mut xnew);
            let rho_new = self.evaluate(&xnew);
            if rho_new < rho {
                rho = rho_new;
                x = xnew;
            } else {
                break;
            }
        }
        (rho, x)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ModalAnalysis
// ─────────────────────────────────────────────────────────────────────────────

/// Result of a modal analysis.
#[derive(Debug, Clone)]
pub struct ModalAnalysisResult {
    /// Natural angular frequencies ωᵢ (rad/s), ascending.
    pub natural_frequencies: Vec<f64>,
    /// Mode shapes φᵢ, each of length n_dof.
    pub mode_shapes: Vec<Vec<f64>>,
    /// Modal masses mᵢ = φᵢᵀ M φᵢ.
    pub modal_masses: Vec<f64>,
    /// Modal stiffnesses kᵢ = φᵢᵀ K φᵢ.
    pub modal_stiffnesses: Vec<f64>,
    /// Modal participation factors Γᵢ.
    pub participation_factors: Vec<f64>,
}

/// Modal analysis driver for FEM structural models.
///
/// Extracts natural frequencies and mode shapes from the generalized eigenvalue
/// problem K·φ = ω²·M·φ and computes derived modal quantities.
pub struct ModalAnalysis {
    /// Row-major n×n stiffness matrix.
    pub k_mat: Vec<f64>,
    /// Row-major n×n mass matrix.
    pub m_mat: Vec<f64>,
    /// Matrix dimension.
    pub n: usize,
    /// Eigensolver configuration.
    pub config: EigensolverConfig,
}

impl ModalAnalysis {
    /// Create a modal analysis problem.
    pub fn new(k_mat: Vec<f64>, m_mat: Vec<f64>, n: usize, config: EigensolverConfig) -> Self {
        Self {
            k_mat,
            m_mat,
            n,
            config,
        }
    }

    /// Run modal analysis and return the full result.
    pub fn run(&self) -> Result<ModalAnalysisResult, String> {
        let gev = GeneralizedEigenvalue::new(
            self.k_mat.clone(),
            self.m_mat.clone(),
            self.n,
            self.config.clone(),
        );
        let (eigenvalues, mode_shapes) = gev.solve_smallest()?;
        let n = self.n;

        let mut natural_frequencies = Vec::with_capacity(eigenvalues.len());
        let mut modal_masses = Vec::with_capacity(eigenvalues.len());
        let mut modal_stiffnesses = Vec::with_capacity(eigenvalues.len());
        let mut participation_factors = Vec::with_capacity(eigenvalues.len());

        // Influence vector (unit vector for base excitation)
        let r: Vec<f64> = vec![1.0; n];

        for (i, phi) in mode_shapes.iter().enumerate() {
            let omega = eigenvalues[i].max(0.0).sqrt();
            natural_frequencies.push(omega);

            let mphi = dense_matvec(&self.m_mat, phi, n);
            let kphi = dense_matvec(&self.k_mat, phi, n);
            let mm: f64 = phi.iter().zip(mphi.iter()).map(|(a, b)| a * b).sum();
            let km: f64 = phi.iter().zip(kphi.iter()).map(|(a, b)| a * b).sum();
            modal_masses.push(mm);
            modal_stiffnesses.push(km);

            // Γᵢ = φᵢᵀ M r / φᵢᵀ M φᵢ
            let mr = dense_matvec(&self.m_mat, &r, n);
            let phimr: f64 = phi.iter().zip(mr.iter()).map(|(a, b)| a * b).sum();
            participation_factors.push(if mm.abs() > 1e-300 { phimr / mm } else { 0.0 });
        }

        Ok(ModalAnalysisResult {
            natural_frequencies,
            mode_shapes,
            modal_masses,
            modal_stiffnesses,
            participation_factors,
        })
    }

    /// Natural frequency in Hz for a given mode index.
    pub fn frequency_hz(omega: f64) -> f64 {
        omega / (2.0 * std::f64::consts::PI)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// FrequencyResponse
// ─────────────────────────────────────────────────────────────────────────────

/// Steady-state frequency-response function (FRF) with hysteretic (structural)
/// damping.
///
/// For a harmonic load `F·e^{iωt}` the complex displacement amplitude is
/// `u(ω) = (K·(1 + 2iζ) − ω²M)⁻¹ F`, i.e. the stiffness carries the loss factor
/// `η = 2ζ`.  For a single degree of freedom this reproduces the viscous
/// resonance amplification `1/(2ζ)`.  `ζ = 0` gives the undamped response
/// `(K − ω²M)⁻¹ F`.  The dynamic stiffness is factorised directly by a dense LU
/// with partial pivoting.
pub struct FrequencyResponse {
    /// Row-major n×n stiffness matrix K.
    pub k_mat: Vec<f64>,
    /// Row-major n×n mass matrix M.
    pub m_mat: Vec<f64>,
    /// Matrix dimension.
    pub n: usize,
    /// Structural damping ratio ζ (0 = undamped); the loss factor is `η = 2ζ`.
    pub damping_ratio: f64,
}

impl FrequencyResponse {
    /// Create a frequency-response solver.
    pub fn new(k_mat: Vec<f64>, m_mat: Vec<f64>, n: usize, damping_ratio: f64) -> Self {
        Self {
            k_mat,
            m_mat,
            n,
            damping_ratio,
        }
    }

    /// Loss factor `η = 2ζ` of the hysteretic damping model.
    pub fn loss_factor(&self) -> f64 {
        2.0 * self.damping_ratio
    }

    /// Complex response `u(ω) = (K(1 + iη) − ω²M)⁻¹ F` as `(real, imaginary)`
    /// parts, with `η = 2ζ`.
    ///
    /// The complex system is solved as the real `2n × 2n` system
    /// `[[K − ω²M, −ηK], [ηK, K − ω²M]]·[u_re; u_im] = [F; 0]`, which reduces to
    /// a single `n × n` solve for an undamped structure.
    ///
    /// # Errors
    ///
    /// Returns an error if `force`, `K` or `M` do not match `n`, if `ω` is not
    /// finite, if the damping ratio is negative or not finite, or if the dynamic
    /// stiffness matrix is singular at this frequency (undamped resonance).
    pub fn response_complex(
        &self,
        omega: f64,
        force: &[f64],
    ) -> Result<(Vec<f64>, Vec<f64>), String> {
        let n = self.n;
        if force.len() != n || self.k_mat.len() != n * n || self.m_mat.len() != n * n {
            return Err(format!(
                "force ({}), K ({}) and M ({}) must match the system size {n}",
                force.len(),
                self.k_mat.len(),
                self.m_mat.len()
            ));
        }
        if !omega.is_finite() {
            return Err("angular frequency must be finite".to_string());
        }
        if !self.damping_ratio.is_finite() || self.damping_ratio < 0.0 {
            return Err("damping ratio must be finite and non-negative".to_string());
        }
        let singular = || "dynamic stiffness matrix is singular at this frequency".to_string();
        let omega_sq = omega * omega;
        // Form dynamic stiffness D = K - ω²M
        let d: Vec<f64> = self
            .k_mat
            .iter()
            .zip(self.m_mat.iter())
            .map(|(k, m)| k - omega_sq * m)
            .collect();
        let eta = self.loss_factor();
        if eta == 0.0 {
            let mut lu = d;
            let perm = lu_factor(&mut lu, n, FRF_PIVOT_TOL).ok_or_else(singular)?;
            return Ok((lu_solve(&lu, &perm, force, n), vec![0.0; n]));
        }
        let big = 2 * n;
        let mut a = vec![0.0_f64; big * big];
        for i in 0..n {
            for j in 0..n {
                let d_ij = d[i * n + j];
                let loss = eta * self.k_mat[i * n + j];
                a[i * big + j] = d_ij;
                a[i * big + n + j] = -loss;
                a[(n + i) * big + j] = loss;
                a[(n + i) * big + n + j] = d_ij;
            }
        }
        let perm = lu_factor(&mut a, big, FRF_PIVOT_TOL).ok_or_else(singular)?;
        let mut rhs = vec![0.0_f64; big];
        rhs[..n].copy_from_slice(force);
        let x = lu_solve(&a, &perm, &rhs, big);
        Ok((x[..n].to_vec(), x[n..].to_vec()))
    }

    /// Response amplitude `|u_i(ω)|` of every DOF for the harmonic load `force`
    /// at angular frequency `omega` (rad/s).
    ///
    /// This is the modulus of the complex response of
    /// [`FrequencyResponse::response_complex`], so the sign (phase) of the
    /// undamped response is not retained; use `response_complex` when the phase
    /// matters.
    ///
    /// Returns `Err` if the dynamic stiffness matrix is singular.
    pub fn response(&self, omega: f64, force: &[f64]) -> Result<Vec<f64>, String> {
        let (re, im) = self.response_complex(omega, force)?;
        Ok(re.iter().zip(im.iter()).map(|(r, i)| r.hypot(*i)).collect())
    }

    /// Compute response over a frequency sweep (rad/s).
    ///
    /// Returns a `Vec<(omega, response_vec)>` for each frequency point.
    pub fn sweep(&self, omegas: &[f64], force: &[f64]) -> Vec<(f64, Result<Vec<f64>, String>)> {
        omegas
            .iter()
            .map(|&om| (om, self.response(om, force)))
            .collect()
    }

    /// Compute the peak response norm over a frequency sweep.
    pub fn peak_response_norm(&self, omegas: &[f64], force: &[f64]) -> f64 {
        self.sweep(omegas, force)
            .into_iter()
            .filter_map(|(_, r)| r.ok())
            .map(|u| norm2(&u))
            .fold(0.0_f64, f64::max)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SpectralNorm
// ─────────────────────────────────────────────────────────────────────────────

/// Computes the spectral norm (largest singular value) of a dense matrix
/// via power iteration on AᵀA.
pub struct SpectralNorm {
    /// Row-major n×n matrix.
    pub a: Vec<f64>,
    /// Matrix dimension.
    pub n: usize,
}

impl SpectralNorm {
    /// Create a spectral norm estimator.
    pub fn new(a: Vec<f64>, n: usize) -> Self {
        assert_eq!(a.len(), n * n);
        Self { a, n }
    }

    /// Estimate the spectral norm ‖A‖₂ (largest singular value).
    ///
    /// Uses power iteration on AᵀA; returns √(dominant eigenvalue of AᵀA).
    pub fn estimate(&self, max_iter: usize, tol: f64) -> f64 {
        let n = self.n;
        // Form AᵀA
        let mut ata = vec![0.0_f64; n * n];
        for i in 0..n {
            for j in 0..n {
                let mut s = 0.0_f64;
                for k in 0..n {
                    s += self.a[k * n + i] * self.a[k * n + j];
                }
                ata[i * n + j] = s;
            }
        }
        let config = EigensolverConfig {
            max_iter,
            tolerance: tol,
            n_modes: 1,
            ..Default::default()
        };
        let pi = PowerIteration::new(ata, n, config);
        let (lam, _) = pi.dominant();
        lam.max(0.0).sqrt()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ModalTruncation
// ─────────────────────────────────────────────────────────────────────────────

/// Reduced-order model via modal truncation.
///
/// Projects the system onto the subspace spanned by the first `n_modes`
/// retained mode shapes, then reconstructs the physical response.
pub struct ModalTruncation {
    /// Retained mode shapes (each of length n_dof).
    pub mode_shapes: Vec<Vec<f64>>,
    /// Corresponding natural frequencies ωᵢ.
    pub natural_frequencies: Vec<f64>,
    /// Modal masses.
    pub modal_masses: Vec<f64>,
    /// Physical DOF count.
    pub n_dof: usize,
}

impl ModalTruncation {
    /// Create a modal truncation model.
    pub fn new(
        mode_shapes: Vec<Vec<f64>>,
        natural_frequencies: Vec<f64>,
        modal_masses: Vec<f64>,
        n_dof: usize,
    ) -> Self {
        Self {
            mode_shapes,
            natural_frequencies,
            modal_masses,
            n_dof,
        }
    }

    /// Project a physical force vector onto modal coordinates.
    ///
    /// Returns modal forces fᵢ = φᵢᵀ · f.
    pub fn modal_force(&self, f: &[f64]) -> Vec<f64> {
        self.mode_shapes.iter().map(|phi| dot(phi, f)).collect()
    }

    /// Compute static response in modal coordinates: qᵢ = fᵢ / kᵢ (kᵢ = ωᵢ² mᵢ).
    pub fn static_modal_response(&self, f: &[f64]) -> Vec<f64> {
        let modal_f = self.modal_force(f);
        modal_f
            .iter()
            .enumerate()
            .map(|(i, &fi)| {
                let ki = self.natural_frequencies[i].powi(2) * self.modal_masses[i];
                if ki.abs() > 1e-300 { fi / ki } else { 0.0 }
            })
            .collect()
    }

    /// Reconstruct physical displacement from modal coordinates q.
    ///
    /// u ≈ Σᵢ qᵢ φᵢ
    pub fn reconstruct(&self, modal_q: &[f64]) -> Vec<f64> {
        let n = self.n_dof;
        let mut u = vec![0.0_f64; n];
        for (i, phi) in self.mode_shapes.iter().enumerate() {
            let qi = if i < modal_q.len() { modal_q[i] } else { 0.0 };
            for (j, &phi_j) in phi.iter().enumerate() {
                u[j] += qi * phi_j;
            }
        }
        u
    }

    /// Compute physical response to static load via modal superposition.
    pub fn static_response(&self, f: &[f64]) -> Vec<f64> {
        let q = self.static_modal_response(f);
        self.reconstruct(&q)
    }

    /// Effective modal mass fraction for mode i: mᵢ_eff / m_total.
    pub fn modal_mass_fraction(&self, m_total: f64) -> Vec<f64> {
        self.modal_masses
            .iter()
            .map(|&mi| mi / m_total.max(1e-300))
            .collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests;
