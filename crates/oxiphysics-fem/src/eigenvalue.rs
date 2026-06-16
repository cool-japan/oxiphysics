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
//! - [`InverseIteration`] — smallest eigenvalue via inverse iteration (LU)
//! - [`LanczosMethod`] — Lanczos tridiagonalization for k eigenvalues
//! - [`GeneralizedEigenvalue`] — K·φ = λM·φ via Cholesky of M
//! - [`ModalAnalysis`] — natural frequencies, mode shapes, modal participation
//! - [`FrequencyResponse`] — steady-state harmonic response
//! - [`RayleighQuotient`] — upper bound on fundamental frequency
//! - [`EigensolverConfig`] — iteration / convergence configuration
//! - [`SpectralNorm`] — largest singular value via power iteration
//! - [`ModalTruncation`] — reduced-order model from retained modes

pub mod lobpcg;
pub use lobpcg::{lobpcg_solve, LobpcgConfig, LobpcgResult, EigensolverError};

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

/// Dense symmetric tridiagonal eigenvalues via implicit QL algorithm (0-based).
/// Returns eigenvalues in ascending order for an n×n symmetric tridiagonal
/// matrix with diagonal `diag` (length n) and sub-diagonal `off` (length n-1).
///
/// Adapted from the standard Numerical Recipes / LAPACK `dstev` approach.
fn tridiagonal_eigenvalues(diag: &[f64], off: &[f64]) -> Vec<f64> {
    let n = diag.len();
    if n == 0 {
        return vec![];
    }
    if n == 1 {
        return vec![diag[0]];
    }

    let mut d = diag.to_vec();
    // e[i] = off-diagonal below row i, e[0] = 0 (unused), e[i] = off[i-1] for i>=1
    let mut e = vec![0.0_f64; n];
    for i in 1..n {
        e[i] = if (i - 1) < off.len() { off[i - 1] } else { 0.0 };
    }

    // Iterate over each unreduced sub-matrix
    for l in 0..n {
        let mut iter = 0usize;
        loop {
            // Find the largest m ∈ [l, n-1) such that e[m+1] is small.
            // Convention: m is the index of the last diagonal in the active block.
            let mut m = l;
            while m < n - 1 {
                let thresh = 1e-15 * (d[m].abs() + d[m + 1].abs());
                if e[m + 1].abs() <= thresh {
                    break;
                }
                m += 1;
            }
            if m == l {
                break;
            } // already diagonal
            if iter >= 30 * n {
                break;
            } // give up
            iter += 1;

            // Wilkinson shift: use the eigenvalue of the bottom 2×2
            // sub-matrix (d[m-1], e[m]; e[m], d[m]) closer to d[m]
            let shift = {
                let b = (d[m - 1] - d[m]) / 2.0;
                let sign_b = if b >= 0.0 { 1.0_f64 } else { -1.0_f64 };
                d[m] - e[m] * e[m] / (b + sign_b * b.hypot(e[m]))
            };

            // Chase the bulge from index l to m
            let mut g = d[l] - shift;
            let mut s = 1.0_f64;
            let mut c = 1.0_f64;
            let mut p = 0.0_f64;

            for i in l..m {
                let f = s * e[i + 1];
                let bb = c * e[i + 1];
                let r = f.hypot(g);
                e[i] = r; // Note: e[l] will be overwritten at the end
                if r.abs() < 1e-300 {
                    d[i] -= p;
                    e[m] = 0.0;
                    break;
                }
                s = f / r;
                c = g / r;
                g = d[i] - p;
                let r2 = (d[i + 1] - g) * s + 2.0 * c * bb;
                p = s * r2;
                d[i] = g + p;
                g = c * r2 - bb;
            }
            d[m] -= p;
            e[m] = g;
            e[l] = 0.0;
        }
    }
    d.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    d
}

/// In-place dense LU factorization (no pivoting) for an n×n row-major matrix.
/// Returns `Err(())` if a zero pivot is encountered.
fn lu_factor(a: &mut [f64], n: usize) -> Result<(), ()> {
    for k in 0..n {
        let pivot = a[k * n + k];
        if pivot.abs() < 1e-300 {
            return Err(());
        }
        for i in (k + 1)..n {
            a[i * n + k] /= pivot;
            for j in (k + 1)..n {
                let m = a[i * n + k];
                a[i * n + j] -= m * a[k * n + j];
            }
        }
    }
    Ok(())
}

/// Solve L·U·x = b in-place.  `lu` contains the factored matrix.
fn lu_solve(lu: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    let mut x = b.to_vec();
    // Forward substitution L·y = b
    for i in 1..n {
        for j in 0..i {
            x[i] -= lu[i * n + j] * x[j];
        }
    }
    // Back substitution U·x = y
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
/// Uses LU factorization of (A − σ·I) and iterates with back-substitution.
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
        lu_factor(&mut lu, n).map_err(|_| "singular shifted matrix".to_string())?;

        // Non-axis-aligned starting vector
        let scale = (n as f64).sqrt().recip();
        let mut v: Vec<f64> = (0..n).map(|i| (i + 1) as f64 * scale).collect();
        normalize(&mut v);
        let mut prev_v = v.clone();

        for _ in 0..self.config.max_iter {
            let w = lu_solve(&lu, &v, n);
            let mut new_v = w;
            normalize(&mut new_v);
            // Check convergence: angle between successive iterates
            let angle_cos = dot(&new_v, &v).abs();
            let converged = (1.0 - angle_cos) < self.config.tolerance;
            prev_v = v;
            v = new_v;
            if converged {
                break;
            }
        }
        let _ = prev_v; // suppress unused warning
        // Rayleigh quotient: λ = vᵀ A v / vᵀ v (vᵀv = 1 since normalized)
        let av = dense_matvec(&self.a, &v, n);
        let eigenvalue = dot(&v, &av);
        Ok((eigenvalue, v))
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// LanczosMethod
// ─────────────────────────────────────────────────────────────────────────────

/// Lanczos tridiagonalization for extracting k eigenvalues of a symmetric matrix.
///
/// Builds a Krylov basis and extracts eigenvalues from the resulting
/// tridiagonal matrix.  Suitable for large sparse problems.
pub struct LanczosMethod {
    /// Row-major n×n symmetric matrix.
    pub a: Vec<f64>,
    /// Matrix dimension.
    pub n: usize,
    /// Solver configuration (`n_modes` determines k).
    pub config: EigensolverConfig,
}

impl LanczosMethod {
    /// Create a new Lanczos solver.
    pub fn new(a: Vec<f64>, n: usize, config: EigensolverConfig) -> Self {
        assert_eq!(a.len(), n * n);
        Self { a, n, config }
    }

    /// Run Lanczos iteration and return up to `config.n_modes` eigenvalues.
    ///
    /// Returns eigenvalues sorted ascending.
    pub fn eigenvalues(&self) -> Vec<f64> {
        let n = self.n;
        let k = self.config.n_modes.min(n);
        // Non-axis-aligned starting vector so all eigenspaces are represented
        let scale = (n as f64).sqrt().recip();
        let mut v0: Vec<f64> = (0..n).map(|i| (i + 1) as f64 * scale).collect();
        normalize(&mut v0);
        let mut q: Vec<Vec<f64>> = Vec::with_capacity(k + 1);
        q.push(v0);

        let mut alpha_vec = Vec::with_capacity(k);
        let mut beta_vec = Vec::with_capacity(k); // beta_vec[j] = β_{j+1}
        let mut steps = 0usize;

        for j in 0..k {
            let w_raw = dense_matvec(&self.a, &q[j], n);
            let a_j: f64 = dot(&w_raw, &q[j]);
            alpha_vec.push(a_j);

            // w = A q_j - a_j q_j - b_j q_{j-1}
            let mut w: Vec<f64> = w_raw
                .iter()
                .enumerate()
                .map(|(idx, &wi)| {
                    let mut val = wi - a_j * q[j][idx];
                    if j > 0 {
                        val -= beta_vec[j - 1] * q[j - 1][idx];
                    }
                    val
                })
                .collect();

            steps += 1;
            let b = norm2(&w);
            beta_vec.push(b);

            if j + 1 < k {
                if b < 1e-14 {
                    // Invariant subspace found — stop early
                    break;
                }
                for x in w.iter_mut() {
                    *x /= b;
                }
                q.push(w);
            }
        }

        let m = steps.min(alpha_vec.len());
        // off-diagonal: beta_vec[0..m-1]
        let off: Vec<f64> = beta_vec[..m.saturating_sub(1)].to_vec();
        let mut eigs = tridiagonal_eigenvalues(&alpha_vec[..m], &off);
        eigs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        eigs.truncate(k);
        eigs
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// GeneralizedEigenvalue
// ─────────────────────────────────────────────────────────────────────────────

/// Generalized eigenvalue solver for K·φ = λM·φ.
///
/// Uses Cholesky factorization of M to convert to a standard eigenvalue problem:
/// `L⁻¹ K L⁻ᵀ ỹ = λ ỹ`, where M = L Lᵀ.
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

    /// Compute the `config.n_modes` smallest generalized eigenvalues.
    ///
    /// Returns `(eigenvalues, eigenvectors)` where each eigenvector has length n.
    pub fn solve_smallest(&self) -> Result<(Vec<f64>, Vec<Vec<f64>>), String> {
        let n = self.n;
        let l = Self::cholesky(&self.m_mat, n)?;

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

        // Use Lanczos on C
        let lanczos = LanczosMethod {
            a: c.clone(),
            n,
            config: self.config.clone(),
        };
        let eigenvalues = lanczos.eigenvalues();

        // Back-transform eigenvectors: φ = L⁻ᵀ ỹ
        // Use power iteration on C to get eigenvectors for the k smallest
        let k = self.config.n_modes.min(n).min(eigenvalues.len());
        let mut eigenvectors = Vec::with_capacity(k);
        for &lam in eigenvalues.iter().take(k) {
            // Inverse iteration on C shifted by lambda
            let inv_cfg = EigensolverConfig {
                shift: lam - 1e-6,
                n_modes: 1,
                ..self.config.clone()
            };
            let inv = InverseIteration {
                a: c.clone(),
                n,
                config: inv_cfg,
            };
            let phi_tilde = match inv.solve() {
                Ok((_l, v)) => v,
                Err(_) => vec![0.0; n],
            };
            let phi = Self::bwd_sub(&l, &phi_tilde, n);
            eigenvectors.push(phi);
        }

        Ok((eigenvalues[..k].to_vec(), eigenvectors))
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

/// Steady-state frequency-response function (FRF).
///
/// Computes the displacement amplitude u(ω) = (K − ω²M)⁻¹ F for a harmonic
/// load F at angular frequency ω, using direct factorization.
pub struct FrequencyResponse {
    /// Row-major n×n stiffness matrix K.
    pub k_mat: Vec<f64>,
    /// Row-major n×n mass matrix M.
    pub m_mat: Vec<f64>,
    /// Matrix dimension.
    pub n: usize,
    /// Structural damping ratio ζ (0 = undamped).
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

    /// Compute real-valued response u(ω) = (K − ω²M)⁻¹ F (undamped or lightly damped).
    ///
    /// Returns `Err` if the dynamic stiffness matrix is singular.
    pub fn response(&self, omega: f64, force: &[f64]) -> Result<Vec<f64>, String> {
        let n = self.n;
        assert_eq!(force.len(), n);
        // Form dynamic stiffness D = K - ω²M
        let mut d = self.k_mat.clone();
        for i in 0..n {
            for j in 0..n {
                d[i * n + j] -= omega * omega * self.m_mat[i * n + j];
            }
        }
        lu_factor(&mut d, n)
            .map_err(|_| "dynamic stiffness matrix is singular at this frequency".to_string())?;
        Ok(lu_solve(&d, force, n))
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
mod tests {
    use super::*;

    // ── Utility helpers ───────────────────────────────────────────────────────

    /// 2×2 symmetric positive-definite matrix: \[\[4,1\\],\[1,3\]]
    fn mat2_spd() -> (Vec<f64>, usize) {
        (vec![4.0, 1.0, 1.0, 3.0], 2)
    }

    /// Identity matrix n×n
    fn identity(n: usize) -> Vec<f64> {
        let mut a = vec![0.0_f64; n * n];
        for i in 0..n {
            a[i * n + i] = 1.0;
        }
        a
    }

    /// Diagonal matrix with given diagonal entries
    fn diag_mat(diag: &[f64]) -> Vec<f64> {
        let n = diag.len();
        let mut a = vec![0.0_f64; n * n];
        for (i, &d) in diag.iter().enumerate() {
            a[i * n + i] = d;
        }
        a
    }

    // ── EigensolverConfig ─────────────────────────────────────────────────────

    #[test]
    fn test_eigensolver_config_default() {
        let cfg = EigensolverConfig::default();
        assert_eq!(cfg.max_iter, 500);
        assert_eq!(cfg.n_modes, 5);
        assert!(cfg.tolerance < 1e-9);
        assert_eq!(cfg.shift, 0.0);
    }

    #[test]
    fn test_eigensolver_config_custom() {
        let cfg = EigensolverConfig {
            max_iter: 100,
            tolerance: 1e-6,
            shift: 2.0,
            n_modes: 3,
        };
        assert_eq!(cfg.n_modes, 3);
        assert_eq!(cfg.shift, 2.0);
    }

    // ── PowerIteration ────────────────────────────────────────────────────────

    #[test]
    fn test_power_iteration_identity() {
        let (a, n) = (identity(3), 3);
        let cfg = EigensolverConfig {
            max_iter: 200,
            tolerance: 1e-10,
            ..Default::default()
        };
        let pi = PowerIteration::new(a, n, cfg);
        let (lam, v) = pi.dominant();
        assert!(
            (lam - 1.0).abs() < 1e-6,
            "eigenvalue of I should be 1, got {lam}"
        );
        assert!((norm2(&v) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_power_iteration_diagonal() {
        let a = diag_mat(&[1.0, 5.0, 3.0]);
        let cfg = EigensolverConfig {
            max_iter: 1000,
            tolerance: 1e-10,
            ..Default::default()
        };
        let pi = PowerIteration::new(a, 3, cfg);
        let (lam, _) = pi.dominant();
        assert!(
            (lam - 5.0).abs() < 1e-5,
            "dominant eigenvalue should be 5, got {lam}"
        );
    }

    #[test]
    fn test_power_iteration_2x2_spd() {
        // [[4,1],[1,3]] has eigenvalues ≈ 4.618, 2.382
        let (a, n) = mat2_spd();
        let cfg = EigensolverConfig {
            max_iter: 500,
            tolerance: 1e-10,
            ..Default::default()
        };
        let pi = PowerIteration::new(a, n, cfg);
        let (lam, v) = pi.dominant();
        assert!(
            lam > 4.0 && lam < 5.0,
            "dominant eigenvalue out of range: {lam}"
        );
        assert!((norm2(&v) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_power_iteration_deflation_count() {
        let a = diag_mat(&[1.0, 4.0, 9.0]);
        let cfg = EigensolverConfig {
            max_iter: 1000,
            tolerance: 1e-10,
            ..Default::default()
        };
        let pi = PowerIteration::new(a, 3, cfg);
        let results = pi.deflated_eigenvalues(2);
        assert_eq!(results.len(), 2);
    }

    #[test]
    fn test_power_iteration_eigvec_unit_norm() {
        let a = diag_mat(&[2.0, 7.0, 3.0]);
        let cfg = EigensolverConfig {
            max_iter: 500,
            tolerance: 1e-10,
            ..Default::default()
        };
        let pi = PowerIteration::new(a, 3, cfg);
        let (_, v) = pi.dominant();
        assert!((norm2(&v) - 1.0).abs() < 1e-9);
    }

    // ── InverseIteration ──────────────────────────────────────────────────────

    #[test]
    fn test_inverse_iteration_identity() {
        let cfg = EigensolverConfig {
            max_iter: 300,
            tolerance: 1e-10,
            shift: 0.0,
            n_modes: 1,
        };
        let inv = InverseIteration::new(identity(3), 3, cfg);
        let (lam, _) = inv.solve().unwrap();
        assert!(
            (lam - 1.0).abs() < 1e-5,
            "smallest eigenvalue of I = 1, got {lam}"
        );
    }

    #[test]
    fn test_inverse_iteration_diagonal() {
        let a = diag_mat(&[3.0, 1.0, 5.0]);
        let cfg = EigensolverConfig {
            max_iter: 500,
            tolerance: 1e-10,
            shift: 0.5,
            n_modes: 1,
        };
        let inv = InverseIteration::new(a, 3, cfg);
        let (lam, _) = inv.solve().unwrap();
        // Nearest eigenvalue to shift=0.5 is 1.0
        assert!((lam - 1.0).abs() < 1e-4, "expected 1.0, got {lam}");
    }

    #[test]
    fn test_inverse_iteration_returns_unit_eigvec() {
        let a = diag_mat(&[2.0, 5.0, 8.0]);
        let cfg = EigensolverConfig {
            max_iter: 300,
            tolerance: 1e-10,
            shift: 1.5,
            n_modes: 1,
        };
        let inv = InverseIteration::new(a, 3, cfg);
        let (_, v) = inv.solve().unwrap();
        assert!((norm2(&v) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_inverse_iteration_2x2() {
        let (a, n) = mat2_spd();
        let cfg = EigensolverConfig {
            max_iter: 500,
            tolerance: 1e-10,
            shift: 2.0,
            n_modes: 1,
        };
        let inv = InverseIteration::new(a, n, cfg);
        let (lam, _) = inv.solve().unwrap();
        // Smallest eigenvalue ≈ 2.382
        assert!(lam > 2.0 && lam < 3.0, "expected ~2.38, got {lam}");
    }

    // ── LanczosMethod ─────────────────────────────────────────────────────────

    #[test]
    fn test_lanczos_identity_eigenvalues() {
        let cfg = EigensolverConfig {
            n_modes: 3,
            ..Default::default()
        };
        let lm = LanczosMethod::new(identity(4), 4, cfg);
        let eigs = lm.eigenvalues();
        for &e in &eigs {
            assert!((e - 1.0).abs() < 0.1, "identity eigenvalue {e} != 1");
        }
    }

    #[test]
    fn test_lanczos_diagonal_sorted() {
        let a = diag_mat(&[5.0, 2.0, 8.0, 1.0]);
        let cfg = EigensolverConfig {
            n_modes: 4,
            max_iter: 300,
            tolerance: 1e-10,
            ..Default::default()
        };
        let lm = LanczosMethod::new(a, 4, cfg);
        let eigs = lm.eigenvalues();
        for w in eigs.windows(2) {
            assert!(w[0] <= w[1] + 1e-10, "eigenvalues not sorted");
        }
    }

    #[test]
    fn test_lanczos_n_modes_limit() {
        let a = diag_mat(&[1.0, 2.0, 3.0, 4.0, 5.0]);
        let cfg = EigensolverConfig {
            n_modes: 3,
            ..Default::default()
        };
        let lm = LanczosMethod::new(a, 5, cfg);
        let eigs = lm.eigenvalues();
        assert!(eigs.len() <= 3);
    }

    #[test]
    fn test_lanczos_2x2_spd() {
        let (a, n) = mat2_spd();
        let cfg = EigensolverConfig {
            n_modes: 2,
            max_iter: 300,
            tolerance: 1e-10,
            ..Default::default()
        };
        let lm = LanczosMethod::new(a, n, cfg);
        let eigs = lm.eigenvalues();
        assert!(!eigs.is_empty());
        // Both eigenvalues should be positive
        for &e in &eigs {
            assert!(e > 0.0, "eigenvalue should be positive, got {e}");
        }
    }

    // ── RayleighQuotient ──────────────────────────────────────────────────────

    #[test]
    fn test_rayleigh_quotient_identity() {
        let k = identity(3);
        let m = identity(3);
        let rq = RayleighQuotient::new(k, m, 3);
        let x = vec![1.0, 0.0, 0.0];
        let r = rq.evaluate(&x);
        assert!((r - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_rayleigh_quotient_diagonal() {
        let k = diag_mat(&[1.0, 4.0, 9.0]);
        let m = identity(3);
        let rq = RayleighQuotient::new(k, m, 3);
        // x aligned with first mode
        let x = vec![1.0, 0.0, 0.0];
        let r = rq.evaluate(&x);
        assert!((r - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_rayleigh_quotient_upper_bound() {
        let k = diag_mat(&[1.0, 4.0, 9.0]);
        let m = identity(3);
        let rq = RayleighQuotient::new(k, m, 3);
        // Any trial vector ρ(x) >= min eigenvalue = 1.0
        let x = vec![0.5, 0.5, 0.5f64.sqrt()];
        let r = rq.evaluate(&x);
        assert!(
            r >= 1.0 - 1e-10,
            "Rayleigh quotient must be >= min eigenvalue"
        );
    }

    #[test]
    fn test_rayleigh_minimise() {
        let k = diag_mat(&[1.0, 4.0, 9.0]);
        let m = identity(3);
        let rq = RayleighQuotient::new(k, m, 3);
        let (rho_min, _) = rq.minimise(100);
        assert!(
            rho_min >= 1.0 - 1e-6,
            "minimised RQ must stay >= min eigenvalue"
        );
        assert!(rho_min <= 9.0 + 1e-6);
    }

    // ── GeneralizedEigenvalue ─────────────────────────────────────────────────

    #[test]
    fn test_generalized_eigenvalue_identity_mass() {
        // K = diag(1,4), M = I → eigenvalues = 1, 4
        let k = diag_mat(&[1.0, 4.0]);
        let m = identity(2);
        let cfg = EigensolverConfig {
            n_modes: 2,
            max_iter: 500,
            tolerance: 1e-10,
            ..Default::default()
        };
        let gev = GeneralizedEigenvalue::new(k, m, 2, cfg);
        let (eigs, _) = gev.solve_smallest().unwrap();
        assert!(!eigs.is_empty());
        // First eigenvalue should be near 1.0
        assert!(
            eigs[0] > 0.5 && eigs[0] < 2.0,
            "expected ~1, got {}",
            eigs[0]
        );
    }

    #[test]
    fn test_generalized_eigenvalue_cholesky_spd() {
        // 2×2 SPD mass matrix
        let m = vec![4.0, 2.0, 2.0, 3.0];
        let l = GeneralizedEigenvalue::cholesky(&m, 2).unwrap();
        // L[0,0] should be sqrt(4) = 2
        assert!((l[0] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_generalized_cholesky_non_spd_fails() {
        let m = vec![-1.0, 0.0, 0.0, 1.0];
        assert!(GeneralizedEigenvalue::cholesky(&m, 2).is_err());
    }

    // ── ModalAnalysis ─────────────────────────────────────────────────────────

    #[test]
    fn test_modal_analysis_frequency_hz() {
        let omega = 2.0 * std::f64::consts::PI; // 1 Hz
        let f = ModalAnalysis::frequency_hz(omega);
        assert!((f - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_modal_analysis_run_counts() {
        let k = diag_mat(&[1.0, 4.0]);
        let m = identity(2);
        let cfg = EigensolverConfig {
            n_modes: 2,
            max_iter: 500,
            tolerance: 1e-10,
            ..Default::default()
        };
        let ma = ModalAnalysis::new(k, m, 2, cfg);
        let result = ma.run().unwrap();
        assert!(result.natural_frequencies.len() <= 2);
        assert_eq!(result.natural_frequencies.len(), result.mode_shapes.len());
        assert_eq!(result.natural_frequencies.len(), result.modal_masses.len());
    }

    #[test]
    fn test_modal_analysis_nonnegative_frequencies() {
        let k = diag_mat(&[1.0, 4.0, 9.0]);
        let m = identity(3);
        let cfg = EigensolverConfig {
            n_modes: 3,
            max_iter: 500,
            tolerance: 1e-10,
            ..Default::default()
        };
        let ma = ModalAnalysis::new(k, m, 3, cfg);
        let result = ma.run().unwrap();
        for &omega in &result.natural_frequencies {
            assert!(omega >= 0.0, "frequency must be non-negative: {omega}");
        }
    }

    // ── FrequencyResponse ─────────────────────────────────────────────────────

    #[test]
    fn test_frequency_response_static_limit() {
        // At ω=0, response = K⁻¹ F  (static solution)
        let k = diag_mat(&[2.0, 8.0]);
        let m = identity(2);
        let frf = FrequencyResponse::new(k, m, 2, 0.0);
        let f = vec![2.0, 8.0];
        let u = frf.response(0.0, &f).unwrap();
        assert!((u[0] - 1.0).abs() < 1e-10);
        assert!((u[1] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_frequency_response_singular_at_resonance() {
        // K = diag(1), M = I → resonance at ω=1
        let k = diag_mat(&[1.0]);
        let m = identity(1);
        let frf = FrequencyResponse::new(k, m, 1, 0.0);
        let result = frf.response(1.0, &[1.0]);
        assert!(result.is_err(), "should fail at exact resonance");
    }

    #[test]
    fn test_frequency_response_sweep_length() {
        let k = diag_mat(&[4.0, 9.0]);
        let m = identity(2);
        let frf = FrequencyResponse::new(k, m, 2, 0.0);
        let omegas: Vec<f64> = (1..=10).map(|i| i as f64 * 0.1).collect();
        let sweep = frf.sweep(&omegas, &[1.0, 1.0]);
        assert_eq!(sweep.len(), 10);
    }

    #[test]
    fn test_frequency_response_peak_norm_positive() {
        let k = diag_mat(&[4.0, 9.0]);
        let m = identity(2);
        let frf = FrequencyResponse::new(k, m, 2, 0.0);
        let omegas: Vec<f64> = (1..=20).map(|i| i as f64 * 0.1).collect();
        let peak = frf.peak_response_norm(&omegas, &[1.0, 1.0]);
        assert!(peak > 0.0);
    }

    // ── SpectralNorm ──────────────────────────────────────────────────────────

    #[test]
    fn test_spectral_norm_identity() {
        let sn = SpectralNorm::new(identity(3), 3);
        let norm = sn.estimate(500, 1e-10);
        assert!(
            (norm - 1.0).abs() < 0.1,
            "spectral norm of I = 1, got {norm}"
        );
    }

    #[test]
    fn test_spectral_norm_diagonal() {
        // diag(1,2,3) — spectral norm = 3
        let a = diag_mat(&[1.0, 2.0, 3.0]);
        let sn = SpectralNorm::new(a, 3);
        let norm = sn.estimate(500, 1e-10);
        assert!(
            (norm - 3.0).abs() < 0.5,
            "spectral norm should be ~3, got {norm}"
        );
    }

    #[test]
    fn test_spectral_norm_nonnegative() {
        let (a, n) = mat2_spd();
        let sn = SpectralNorm::new(a, n);
        let norm = sn.estimate(200, 1e-8);
        assert!(norm >= 0.0);
    }

    // ── ModalTruncation ───────────────────────────────────────────────────────

    #[test]
    fn test_modal_truncation_reconstruct_zero() {
        let modes = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let freqs = vec![1.0, 2.0];
        let masses = vec![1.0, 1.0];
        let mt = ModalTruncation::new(modes, freqs, masses, 2);
        let u = mt.reconstruct(&[0.0, 0.0]);
        assert_eq!(u, vec![0.0, 0.0]);
    }

    #[test]
    fn test_modal_truncation_reconstruct_first_mode() {
        let modes = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let freqs = vec![1.0, 2.0];
        let masses = vec![1.0, 1.0];
        let mt = ModalTruncation::new(modes, freqs, masses, 2);
        let u = mt.reconstruct(&[3.0, 0.0]);
        assert!((u[0] - 3.0).abs() < 1e-12);
        assert!(u[1].abs() < 1e-12);
    }

    #[test]
    fn test_modal_truncation_static_response() {
        // K = diag(1,4), modes = unit, frequencies = 1,2, masses = 1,1
        let modes = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let freqs = vec![1.0, 2.0];
        let masses = vec![1.0, 1.0];
        let mt = ModalTruncation::new(modes, freqs, masses, 2);
        // Force [1,0] → modal force = [1,0], modal k = [1,4]
        // q = [1/1, 0/4] = [1, 0] → u = [1, 0]
        let u = mt.static_response(&[1.0, 0.0]);
        assert!((u[0] - 1.0).abs() < 1e-10);
        assert!(u[1].abs() < 1e-10);
    }

    #[test]
    fn test_modal_truncation_mass_fraction() {
        let modes = vec![vec![1.0, 0.0]];
        let freqs = vec![1.0];
        let masses = vec![3.0];
        let mt = ModalTruncation::new(modes, freqs, masses, 2);
        let frac = mt.modal_mass_fraction(10.0);
        assert!((frac[0] - 0.3).abs() < 1e-10);
    }

    #[test]
    fn test_modal_truncation_modal_force() {
        let modes = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let freqs = vec![1.0, 2.0];
        let masses = vec![1.0, 1.0];
        let mt = ModalTruncation::new(modes, freqs, masses, 2);
        let mf = mt.modal_force(&[5.0, 7.0]);
        assert!((mf[0] - 5.0).abs() < 1e-12);
        assert!((mf[1] - 7.0).abs() < 1e-12);
    }

    // ── Utility functions ─────────────────────────────────────────────────────

    #[test]
    fn test_dot_product() {
        assert!((dot(&[1.0, 2.0, 3.0], &[4.0, 5.0, 6.0]) - 32.0).abs() < 1e-12);
    }

    #[test]
    fn test_norm2() {
        assert!((norm2(&[3.0, 4.0]) - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_normalize() {
        let mut v = vec![3.0, 4.0];
        let n = normalize(&mut v);
        assert!((n - 5.0).abs() < 1e-12);
        assert!((norm2(&v) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_dense_matvec_identity() {
        let a = identity(3);
        let x = vec![1.0, 2.0, 3.0];
        let y = dense_matvec(&a, &x, 3);
        assert_eq!(y, x);
    }

    #[test]
    fn test_lu_factor_solve() {
        // [[2,1],[4,3]] x = [1,1] → x = [1,-1]
        let mut a = vec![2.0, 1.0, 4.0, 3.0];
        lu_factor(&mut a, 2).unwrap();
        let x = lu_solve(&a, &[1.0, 1.0], 2);
        assert!((x[0] - 1.0).abs() < 1e-10, "x[0]={}", x[0]);
        assert!((x[1] + 1.0).abs() < 1e-10, "x[1]={}", x[1]);
    }

    #[test]
    fn test_lu_factor_identity() {
        let mut a = identity(3);
        let res = lu_factor(&mut a, 3);
        assert!(res.is_ok());
    }

    #[test]
    fn test_tridiagonal_eigenvalues_sorted() {
        let d = vec![2.0, 4.0, 6.0];
        let e = vec![0.5, 0.5];
        let eigs = tridiagonal_eigenvalues(&d, &e);
        for w in eigs.windows(2) {
            assert!(w[0] <= w[1]);
        }
    }

    #[test]
    fn test_tridiagonal_eigenvalues_diagonal_case() {
        // Zero off-diagonal → eigenvalues = diagonal
        let d = vec![1.0, 3.0, 2.0];
        let e = vec![0.0, 0.0];
        let eigs = tridiagonal_eigenvalues(&d, &e);
        assert_eq!(eigs.len(), 3);
        assert!((eigs[0] - 1.0).abs() < 1e-8);
        assert!((eigs[1] - 2.0).abs() < 1e-8);
        assert!((eigs[2] - 3.0).abs() < 1e-8);
    }
}
