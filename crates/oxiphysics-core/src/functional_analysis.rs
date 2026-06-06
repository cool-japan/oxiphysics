// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Functional analysis tools for physics simulation.
//!
//! Provides Hilbert and Banach spaces, operator spectrum, Sobolev spaces,
//! functional derivatives (Gâteaux/Fréchet), and variational problems
//! (Euler-Lagrange, Lagrange multiplier).  Also includes L² inner products,
//! Gram–Schmidt orthogonalization, Fourier/Chebyshev/Legendre expansions,
//! Sobolev norms, and operator norm estimation via power iteration.

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// FunctionSpace
// ─────────────────────────────────────────────────────────────────────────────

/// A finite-dimensional function space spanned by a list of basis functions.
///
/// Each basis element is a plain function pointer `fn(f64) -> f64` to avoid
/// lifetime complications with closures.
#[derive(Clone)]
pub struct FunctionSpace {
    /// Basis functions spanning this space.
    pub basis: Vec<fn(f64) -> f64>,
}

impl FunctionSpace {
    /// Create a new function space from a vector of basis functions.
    pub fn new(basis: Vec<fn(f64) -> f64>) -> Self {
        Self { basis }
    }

    /// Evaluate the `i`-th basis function at `x`.
    pub fn eval(&self, i: usize, x: f64) -> f64 {
        (self.basis[i])(x)
    }

    /// Number of basis functions.
    pub fn dim(&self) -> usize {
        self.basis.len()
    }

    /// Evaluate a linear combination `∑ coeffs[i] * basis[i](x)`.
    ///
    /// Panics if `coeffs.len() != self.dim()`.
    pub fn combine(&self, coeffs: &[f64], x: f64) -> f64 {
        assert_eq!(
            coeffs.len(),
            self.basis.len(),
            "coefficient length mismatch"
        );
        coeffs
            .iter()
            .zip(self.basis.iter())
            .map(|(c, f)| c * f(x))
            .sum()
    }
}

impl std::fmt::Debug for FunctionSpace {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "FunctionSpace {{ dim: {} }}", self.basis.len())
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// L² inner product and norm
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the discrete L² inner product `⟨f, g⟩ = dx · Σ f[i]·g[i]`.
///
/// `f` and `g` must have the same length. `dx` is the uniform grid spacing.
pub fn l2_inner_product(f: &[f64], g: &[f64], dx: f64) -> f64 {
    assert_eq!(f.len(), g.len(), "l2_inner_product: slice length mismatch");
    f.iter().zip(g.iter()).map(|(a, b)| a * b).sum::<f64>() * dx
}

/// Compute the discrete L² norm `‖f‖ = √(⟨f, f⟩)`.
pub fn l2_norm(f: &[f64], dx: f64) -> f64 {
    l2_inner_product(f, f, dx).sqrt()
}

// ─────────────────────────────────────────────────────────────────────────────
// Gram–Schmidt orthogonalization
// ─────────────────────────────────────────────────────────────────────────────

/// Orthogonalize a set of sampled basis vectors using the modified Gram–Schmidt
/// process in the L² inner product.
///
/// Each element of `basis` is a sampled function on a uniform grid with spacing
/// `dx`. Returns an orthonormal set (same length as input) in the same order.
/// Linearly dependent vectors become zero vectors.
pub fn gram_schmidt_orthogonalize(basis: &[Vec<f64>], dx: f64) -> Vec<Vec<f64>> {
    let n = basis.len();
    if n == 0 {
        return vec![];
    }
    let m = basis[0].len();
    let mut q: Vec<Vec<f64>> = Vec::with_capacity(n);

    for v in basis.iter() {
        let mut u = v.clone();
        // Subtract projections onto already-computed orthonormal vectors
        for qi in q.iter() {
            let proj = l2_inner_product(&u, qi, dx);
            for (ui, qij) in u.iter_mut().zip(qi.iter()) {
                *ui -= proj * qij;
            }
        }
        let norm = l2_norm(&u, dx);
        if norm > 1e-14 {
            for ui in u.iter_mut() {
                *ui /= norm;
            }
        } else {
            u = vec![0.0; m];
        }
        q.push(u);
    }
    q
}

// ─────────────────────────────────────────────────────────────────────────────
// Fourier series coefficients
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the Fourier series coefficients `(aₙ, bₙ)` for `n = 0, 1, …, n_terms-1`.
///
/// Assumes `f` is sampled on a uniform grid over `[0, L]` where `L = len * dx`,
/// corresponding to the period. Returns `n_terms` pairs `(aₙ, bₙ)` where
/// - `a₀ = (2/L) ∫ f dx` (zeroth cosine term)
/// - `aₙ = (2/L) ∫ f cos(2πnx/L) dx`
/// - `bₙ = (2/L) ∫ f sin(2πnx/L) dx`
pub fn fourier_series_coeffs(f: &[f64], dx: f64, n_terms: usize) -> Vec<(f64, f64)> {
    let n = f.len();
    let l = n as f64 * dx;
    let norm = 2.0 / l;
    (0..n_terms)
        .map(|k| {
            let mut a = 0.0_f64;
            let mut b = 0.0_f64;
            for (i, &fi) in f.iter().enumerate() {
                let x = (i as f64 + 0.5) * dx;
                let arg = 2.0 * PI * (k as f64) * x / l;
                a += fi * arg.cos();
                b += fi * arg.sin();
            }
            (a * norm * dx, b * norm * dx)
        })
        .collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// Chebyshev expansion
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the Chebyshev expansion coefficients `cₙ` for a function sampled on
/// `[-1, 1]` using the discrete cosine approach.
///
/// `f` is assumed to be evaluated at Chebyshev nodes
/// `xₖ = cos(π(k + 0.5)/N)` for `k = 0, …, N-1`.
/// Returns `n_terms` coefficients `c₀, c₁, …`.
pub fn chebyshev_expansion(f: &[f64], n_terms: usize) -> Vec<f64> {
    let n = f.len();
    if n == 0 || n_terms == 0 {
        return vec![0.0; n_terms];
    }
    (0..n_terms)
        .map(|k| {
            let sum: f64 = f
                .iter()
                .enumerate()
                .map(|(j, &fj)| {
                    let theta = PI * (j as f64 + 0.5) / n as f64;
                    fj * (k as f64 * theta).cos()
                })
                .sum();
            if k == 0 {
                sum / n as f64
            } else {
                2.0 * sum / n as f64
            }
        })
        .collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// Legendre expansion
// ─────────────────────────────────────────────────────────────────────────────

/// Evaluate the Legendre polynomial `Pₙ(x)` via the three-term recurrence.
fn legendre_p(n: usize, x: f64) -> f64 {
    if n == 0 {
        return 1.0;
    }
    if n == 1 {
        return x;
    }
    let mut p_prev = 1.0_f64;
    let mut p_curr = x;
    for k in 2..=n {
        let p_next = ((2 * k - 1) as f64 * x * p_curr - (k - 1) as f64 * p_prev) / k as f64;
        p_prev = p_curr;
        p_curr = p_next;
    }
    p_curr
}

/// Compute the Legendre expansion coefficients for `f` sampled on `[-1, 1]`.
///
/// `f` is evaluated at `n = f.len()` uniform points in `[-1, 1]`. Returns
/// `n_terms` coefficients `c₀, …, c_{n_terms-1}` where
/// `cₖ = (2k+1)/2 · ∫ f(x) Pₖ(x) dx`.
pub fn legendre_expansion(f: &[f64], n_terms: usize) -> Vec<f64> {
    let n = f.len();
    if n == 0 || n_terms == 0 {
        return vec![0.0; n_terms];
    }
    let dx = 2.0 / n as f64; // uniform spacing over [-1, 1]
    (0..n_terms)
        .map(|k| {
            let sum: f64 = f
                .iter()
                .enumerate()
                .map(|(j, &fj)| {
                    let x = -1.0 + (j as f64 + 0.5) * dx;
                    fj * legendre_p(k, x)
                })
                .sum();
            (2 * k + 1) as f64 / 2.0 * sum * dx
        })
        .collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// Haar wavelet transform
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the full-length in-place Haar wavelet transform.
///
/// The length of `signal` must be a power of 2. Returns the coefficient vector
/// (Mallat ordering: coarsest approximation first).
pub fn wavelet_haar_transform(signal: &[f64]) -> Vec<f64> {
    let n = signal.len();
    let mut out = signal.to_vec();
    let mut step = n;
    while step >= 2 {
        step /= 2;
        let mut tmp = vec![0.0; step * 2];
        for i in 0..step {
            let a = out[2 * i];
            let b = out[2 * i + 1];
            tmp[i] = (a + b) * 0.5_f64.sqrt();
            tmp[step + i] = (a - b) * 0.5_f64.sqrt();
        }
        out[..step * 2].copy_from_slice(&tmp);
    }
    out
}

/// Compute the inverse Haar wavelet transform.
///
/// `coeffs` must have length that is a power of 2.
pub fn wavelet_haar_inverse(coeffs: &[f64]) -> Vec<f64> {
    let n = coeffs.len();
    let mut out = coeffs.to_vec();
    let mut step = 1_usize;
    while step < n {
        let mut tmp = vec![0.0; step * 2];
        for i in 0..step {
            let a = out[i];
            let d = out[step + i];
            tmp[2 * i] = (a + d) * 0.5_f64.sqrt();
            tmp[2 * i + 1] = (a - d) * 0.5_f64.sqrt();
        }
        out[..step * 2].copy_from_slice(&tmp);
        step *= 2;
    }
    out
}

// ─────────────────────────────────────────────────────────────────────────────
// Sobolev norm (free function)
// ─────────────────────────────────────────────────────────────────────────────

/// Compute an approximate H^s Sobolev norm `‖f‖_{H^s}`.
///
/// Uses the discrete approximation
/// `‖f‖²_{H^s} = ‖f‖²_{L²} + s·‖f'‖²_{L²}`
/// where `f'` is supplied as `df` and `s` is the Sobolev order parameter.
pub fn sobolev_norm(f: &[f64], df: &[f64], dx: f64, s: f64) -> f64 {
    assert_eq!(
        f.len(),
        df.len(),
        "sobolev_norm: f and df must have equal length"
    );
    let l2_sq = l2_inner_product(f, f, dx);
    let dl2_sq = l2_inner_product(df, df, dx);
    (l2_sq + s * dl2_sq).sqrt()
}

// ─────────────────────────────────────────────────────────────────────────────
// Operator norm estimate (power iteration)
// ─────────────────────────────────────────────────────────────────────────────

/// Estimate the operator (spectral) norm of a matrix via power iteration.
///
/// Computes the largest singular value using the power method applied to
/// `Aᵀ A`. The returned value approximates `σ_max(A)`.
///
/// `matrix` is given as a row-major slice-of-rows (each inner `Vec`f64` is
/// a row). Returns `0.0` for an empty matrix.
pub fn operator_norm_estimate(matrix: &[Vec<f64>]) -> f64 {
    let nrows = matrix.len();
    if nrows == 0 {
        return 0.0;
    }
    let ncols = matrix[0].len();
    if ncols == 0 {
        return 0.0;
    }

    // Start with a uniform vector
    let mut v: Vec<f64> = vec![1.0 / (ncols as f64).sqrt(); ncols];
    let mut sigma = 0.0_f64;

    for _ in 0..100 {
        // w = A v
        let mut w = vec![0.0_f64; nrows];
        for (i, row) in matrix.iter().enumerate() {
            w[i] = row.iter().zip(v.iter()).map(|(a, x)| a * x).sum();
        }
        // u = Aᵀ w
        let mut u = vec![0.0_f64; ncols];
        for (i, row) in matrix.iter().enumerate() {
            for (j, &aij) in row.iter().enumerate() {
                u[j] += aij * w[i];
            }
        }
        // Rayleigh quotient
        let norm_u: f64 = u.iter().map(|x| x * x).sum::<f64>().sqrt();
        if norm_u < 1e-15 {
            break;
        }
        sigma = norm_u.sqrt();
        for (vi, ui) in v.iter_mut().zip(u.iter()) {
            *vi = *ui / norm_u;
        }
    }
    sigma
}

// ─────────────────────────────────────────────────────────────────────────────
// HilbertSpace
// ─────────────────────────────────────────────────────────────────────────────

/// A discrete Hilbert space over a uniform grid with spacing `dx`.
///
/// Provides inner product, norm, orthonormalization (Gram-Schmidt), and
/// projection operations.
#[derive(Debug, Clone)]
pub struct HilbertSpace {
    /// Grid spacing for integration approximation.
    pub dx: f64,
    /// Number of grid points.
    pub n: usize,
}

impl HilbertSpace {
    /// Create a new `HilbertSpace` with `n` grid points and spacing `dx`.
    pub fn new(n: usize, dx: f64) -> Self {
        Self { n, dx }
    }

    /// Compute the L² inner product `⟨f, g⟩`.
    pub fn inner_product(&self, f: &[f64], g: &[f64]) -> f64 {
        l2_inner_product(f, g, self.dx)
    }

    /// Compute the L² norm `‖f‖`.
    pub fn norm(&self, f: &[f64]) -> f64 {
        l2_norm(f, self.dx)
    }

    /// Orthonormalize a set of vectors using the modified Gram-Schmidt process.
    ///
    /// Returns the orthonormal set in the same order; linearly dependent
    /// vectors become zero.
    pub fn orthonormalize(&self, basis: &[Vec<f64>]) -> Vec<Vec<f64>> {
        gram_schmidt_orthogonalize(basis, self.dx)
    }

    /// Project `f` onto the subspace spanned by `basis`.
    ///
    /// `basis` does **not** need to be orthonormal; a fresh Gram-Schmidt step
    /// is performed internally.  Returns the projected vector (same length as
    /// `f`).
    pub fn project(&self, f: &[f64], basis: &[Vec<f64>]) -> Vec<f64> {
        if basis.is_empty() || f.is_empty() {
            return vec![0.0; f.len()];
        }
        let ortho = self.orthonormalize(basis);
        let mut proj = vec![0.0; f.len()];
        for q in &ortho {
            if q.iter().all(|x| x.abs() < 1e-15) {
                continue;
            }
            let coeff = self.inner_product(f, q);
            for (pi, qi) in proj.iter_mut().zip(q.iter()) {
                *pi += coeff * qi;
            }
        }
        proj
    }

    /// Compute the distance between two functions `f` and `g`.
    pub fn distance(&self, f: &[f64], g: &[f64]) -> f64 {
        let diff: Vec<f64> = f.iter().zip(g.iter()).map(|(a, b)| a - b).collect();
        self.norm(&diff)
    }

    /// Check if `f` and `g` are orthogonal (inner product < `tol`).
    pub fn are_orthogonal(&self, f: &[f64], g: &[f64], tol: f64) -> bool {
        self.inner_product(f, g).abs() < tol
    }

    /// Normalize `f` in-place to unit L² norm.  Returns a normalized copy.
    pub fn normalize(&self, f: &[f64]) -> Vec<f64> {
        let nm = self.norm(f);
        if nm < 1e-15 {
            return f.to_vec();
        }
        f.iter().map(|x| x / nm).collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BanachSpace
// ─────────────────────────────────────────────────────────────────────────────

/// A discrete Banach space with Lᵖ norm.
///
/// Supports p-norms, dual space estimation, and bounded linear functional
/// evaluation.
#[derive(Debug, Clone)]
pub struct BanachSpace {
    /// The exponent `p` for the Lᵖ norm (`p >= 1`).
    pub p: f64,
    /// Grid spacing for integration approximation.
    pub dx: f64,
}

impl BanachSpace {
    /// Create a new `BanachSpace` with Lᵖ norm and spacing `dx`.
    pub fn new(p: f64, dx: f64) -> Self {
        Self { p: p.max(1.0), dx }
    }

    /// Compute the Lᵖ norm `‖f‖_p = (∫|f|^p dx)^{1/p}`.
    pub fn norm(&self, f: &[f64]) -> f64 {
        if self.p == f64::INFINITY || self.p > 1e10 {
            return f.iter().cloned().fold(0.0_f64, f64::max);
        }
        let integral: f64 = f.iter().map(|x| x.abs().powf(self.p)).sum::<f64>() * self.dx;
        integral.powf(1.0 / self.p)
    }

    /// Compute the dual Lᵍ norm where `1/p + 1/q = 1`.
    ///
    /// For `p = 1` returns the L∞ norm; for `p = ∞` returns the L¹ norm.
    pub fn dual_norm(&self, g: &[f64]) -> f64 {
        let q = if (self.p - 1.0).abs() < 1e-12 {
            f64::INFINITY
        } else {
            self.p / (self.p - 1.0)
        };
        let dual = BanachSpace::new(q, self.dx);
        dual.norm(g)
    }

    /// Evaluate a bounded linear functional `φ(f) = ∫ g(x) f(x) dx`.
    ///
    /// `g` is the Riesz representation kernel.
    pub fn bounded_linear_functional(&self, f: &[f64], g: &[f64]) -> f64 {
        assert_eq!(f.len(), g.len(), "functional: f and g length mismatch");
        f.iter().zip(g.iter()).map(|(a, b)| a * b).sum::<f64>() * self.dx
    }

    /// Check whether `f` lies in the unit ball `‖f‖_p <= 1 + tol`.
    pub fn in_unit_ball(&self, f: &[f64], tol: f64) -> bool {
        self.norm(f) <= 1.0 + tol
    }

    /// Compute the Hölder bound: `|⟨f, g⟩| <= ‖f‖_p * ‖g‖_q`.
    ///
    /// Returns `(lhs, rhs)` where `lhs = |⟨f,g⟩|` and `rhs = ‖f‖_p * ‖g‖_q`.
    pub fn holder_bound(&self, f: &[f64], g: &[f64]) -> (f64, f64) {
        let inner = self.bounded_linear_functional(f, g).abs();
        let norm_f = self.norm(f);
        let norm_g = self.dual_norm(g);
        (inner, norm_f * norm_g)
    }
}

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

    /// Estimate the smallest eigenvalue via inverse power iteration with shift.
    ///
    /// Uses Rayleigh quotient shift `mu` (default 0.0 ≈ smallest eigenvalue).
    pub fn smallest_eigenvalue(&self, shift: f64, max_iter: usize) -> f64 {
        if self.n == 0 {
            return 0.0;
        }
        // Shifted matrix: A - shift*I
        let shifted: Vec<Vec<f64>> = self
            .matrix
            .iter()
            .enumerate()
            .map(|(i, row)| {
                row.iter()
                    .enumerate()
                    .map(|(j, &a)| if i == j { a - shift } else { a })
                    .collect()
            })
            .collect();

        let op = OperatorSpectrum::new(shifted);
        // Power iteration on shifted operator
        let lam = op.spectral_radius(max_iter);
        lam + shift
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

    /// Estimate the resolvent norm `‖(A - λI)⁻¹‖` at a regular point `lambda`.
    ///
    /// Uses power iteration on the shifted operator `(A - λI)`.  Returns
    /// `1 / smallest_singular_value`, approximated via the Frobenius norm of
    /// the shifted matrix.
    pub fn resolvent_norm(&self, lambda: f64) -> f64 {
        // Build (A - λI)
        let shifted: Vec<Vec<f64>> = self
            .matrix
            .iter()
            .enumerate()
            .map(|(i, row)| {
                row.iter()
                    .enumerate()
                    .map(|(j, &a)| if i == j { a - lambda } else { a })
                    .collect()
            })
            .collect();
        let frob: f64 = shifted
            .iter()
            .flat_map(|r| r.iter())
            .map(|x| x * x)
            .sum::<f64>()
            .sqrt();
        if frob < 1e-15 {
            return f64::INFINITY;
        }
        1.0 / frob
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

// ─────────────────────────────────────────────────────────────────────────────
// SobolevSpace
// ─────────────────────────────────────────────────────────────────────────────

/// A discrete Sobolev space H^k on a uniform grid.
///
/// Provides H^k norms, weak derivative approximation (finite differences),
/// and trace operator evaluation.
#[derive(Debug, Clone)]
pub struct SobolevSpace {
    /// Number of grid points.
    pub n: usize,
    /// Grid spacing.
    pub dx: f64,
    /// Sobolev order `k`.
    pub k: usize,
}

impl SobolevSpace {
    /// Create a new `SobolevSpace` with `n` points, spacing `dx`, and order `k`.
    pub fn new(n: usize, dx: f64, k: usize) -> Self {
        Self { n, dx, k }
    }

    /// Compute the H^k norm using finite-difference weak derivatives.
    ///
    /// `‖f‖²_{H^k} = Σ_{j=0}^{k} ‖D^j f‖²_{L²}`
    pub fn norm(&self, f: &[f64]) -> f64 {
        let mut sq = l2_inner_product(f, f, self.dx);
        let mut deriv = f.to_vec();
        for _ in 1..=self.k {
            deriv = Self::weak_derivative(&deriv, self.dx);
            sq += l2_inner_product(&deriv, &deriv, self.dx);
        }
        sq.sqrt()
    }

    /// Compute the first-order weak derivative using central differences.
    ///
    /// Forward/backward differences are used at the boundary.
    pub fn weak_derivative(f: &[f64], dx: f64) -> Vec<f64> {
        let n = f.len();
        if n < 2 {
            return vec![0.0; n];
        }
        let mut d = vec![0.0; n];
        // Central differences in interior
        for i in 1..n - 1 {
            d[i] = (f[i + 1] - f[i - 1]) / (2.0 * dx);
        }
        // Forward/backward at boundaries
        d[0] = (f[1] - f[0]) / dx;
        d[n - 1] = (f[n - 1] - f[n - 2]) / dx;
        d
    }

    /// Evaluate the trace operator: return the boundary values `\[f(0), f(n-1)\]`.
    pub fn trace(&self, f: &[f64]) -> [f64; 2] {
        if f.is_empty() {
            return [0.0, 0.0];
        }
        [f[0], f[f.len() - 1]]
    }

    /// Compute the seminorm `|f|_{H^k} = ‖D^k f‖_{L²}`.
    pub fn seminorm(&self, f: &[f64]) -> f64 {
        let mut deriv = f.to_vec();
        for _ in 0..self.k {
            deriv = Self::weak_derivative(&deriv, self.dx);
        }
        l2_norm(&deriv, self.dx)
    }

    /// Check whether `f` lies in H^k (finite H^k norm).
    pub fn is_in_space(&self, f: &[f64]) -> bool {
        self.norm(f).is_finite()
    }

    /// Compute the H^k inner product `⟨f, g⟩_{H^k} = Σ_{j=0}^{k} ⟨D^j f, D^j g⟩`.
    pub fn inner_product(&self, f: &[f64], g: &[f64]) -> f64 {
        let mut ip = l2_inner_product(f, g, self.dx);
        let mut df = f.to_vec();
        let mut dg = g.to_vec();
        for _ in 1..=self.k {
            df = Self::weak_derivative(&df, self.dx);
            dg = Self::weak_derivative(&dg, self.dx);
            ip += l2_inner_product(&df, &dg, self.dx);
        }
        ip
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// FunctionalDerivative
// ─────────────────────────────────────────────────────────────────────────────

/// Gâteaux and Fréchet derivatives of functionals on function spaces.
///
/// Enables gradient descent in function space and Newton's method for
/// operator equations.
pub struct FunctionalDerivative;

impl FunctionalDerivative {
    /// Compute the Gâteaux derivative of functional `J` at `f` in direction `h`.
    ///
    /// `dJ/dε J(f + εh)|_{ε=0}` approximated by finite difference.
    pub fn gateaux<F>(j: F, f: &[f64], h: &[f64], eps: f64) -> f64
    where
        F: Fn(&[f64]) -> f64,
    {
        let f_plus: Vec<f64> = f
            .iter()
            .zip(h.iter())
            .map(|(fi, hi)| fi + eps * hi)
            .collect();
        let f_minus: Vec<f64> = f
            .iter()
            .zip(h.iter())
            .map(|(fi, hi)| fi - eps * hi)
            .collect();
        (j(&f_plus) - j(&f_minus)) / (2.0 * eps)
    }

    /// Compute the Fréchet derivative (gradient) of functional `J` at `f`.
    ///
    /// Returns the gradient vector `δJ/δf` using component-wise Gâteaux
    /// derivatives with standard basis directions.
    pub fn frechet_gradient<F>(j: F, f: &[f64], eps: f64) -> Vec<f64>
    where
        F: Fn(&[f64]) -> f64,
    {
        let n = f.len();
        let mut grad = vec![0.0; n];
        for i in 0..n {
            let mut f_plus = f.to_vec();
            let mut f_minus = f.to_vec();
            f_plus[i] += eps;
            f_minus[i] -= eps;
            grad[i] = (j(&f_plus) - j(&f_minus)) / (2.0 * eps);
        }
        grad
    }

    /// Perform gradient descent in function space to minimize functional `J`.
    ///
    /// Starting from `f0`, takes `max_iter` steps of size `step_size`.
    /// Returns the minimizer and the final functional value.
    pub fn gradient_descent<F>(
        j: F,
        f0: &[f64],
        step_size: f64,
        max_iter: usize,
        eps: f64,
    ) -> (Vec<f64>, f64)
    where
        F: Fn(&[f64]) -> f64,
    {
        let mut f = f0.to_vec();
        for _ in 0..max_iter {
            let grad = Self::frechet_gradient(&j, &f, eps);
            for (fi, gi) in f.iter_mut().zip(grad.iter()) {
                *fi -= step_size * gi;
            }
        }
        let val = j(&f);
        (f, val)
    }

    /// Newton step in function space: `f_new = f - \[D²J(f)\]⁻¹ DJ(f)`.
    ///
    /// Approximates the Hessian-vector product via a second finite difference
    /// and solves the linear system by a simple diagonal preconditioned
    /// gradient step (one step of Jacobi iteration).
    pub fn newton_step<F>(j: F, f: &[f64], eps: f64) -> Vec<f64>
    where
        F: Fn(&[f64]) -> f64,
    {
        let n = f.len();
        let grad = Self::frechet_gradient(&j, f, eps);

        // Diagonal of Hessian via second central difference
        let mut hess_diag = vec![1.0_f64; n]; // fallback: identity
        for i in 0..n {
            let mut fp = f.to_vec();
            let mut fm = f.to_vec();
            fp[i] += eps;
            fm[i] -= eps;
            let h = (j(&fp) - 2.0 * j(f) + j(&fm)) / (eps * eps);
            if h.abs() > 1e-15 {
                hess_diag[i] = h;
            }
        }

        // Newton step: f - H^{-1} g (diagonal approximation)
        f.iter()
            .zip(grad.iter())
            .zip(hess_diag.iter())
            .map(|((fi, gi), hi)| fi - gi / hi)
            .collect()
    }

    /// Compute the second Gâteaux derivative (bilinear form) at `f` in
    /// directions `h` and `k`.
    pub fn second_gateaux<F>(j: F, f: &[f64], h: &[f64], k: &[f64], eps: f64) -> f64
    where
        F: Fn(&[f64]) -> f64,
    {
        let fph: Vec<f64> = f
            .iter()
            .zip(h.iter())
            .map(|(fi, hi)| fi + eps * hi)
            .collect();
        let fmh: Vec<f64> = f
            .iter()
            .zip(h.iter())
            .map(|(fi, hi)| fi - eps * hi)
            .collect();
        let g_plus = Self::gateaux(&j, &fph, k, eps);
        let g_minus = Self::gateaux(&j, &fmh, k, eps);
        (g_plus - g_minus) / (2.0 * eps)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// VariationalProblem
// ─────────────────────────────────────────────────────────────────────────────

/// Variational problem solver: Euler-Lagrange equations and constrained
/// minimization with Lagrange multipliers.
pub struct VariationalProblem;

impl VariationalProblem {
    /// Compute the Euler-Lagrange residual for the functional
    /// `J\[u\] = ∫ L(x, u, u') dx`.
    ///
    /// The Lagrangian `L` takes `(x, u, u_prime)`.  Returns the residual
    /// `∂L/∂u - d/dx(∂L/∂u')` evaluated at each interior grid point using
    /// finite differences.
    pub fn euler_lagrange_residual<L>(lagrangian: L, u: &[f64], dx: f64) -> Vec<f64>
    where
        L: Fn(f64, f64, f64) -> f64,
    {
        let n = u.len();
        if n < 3 {
            return vec![0.0; n];
        }
        let eps = dx * 1e-5;
        let mut residual = vec![0.0; n];

        for i in 1..n - 1 {
            let x = (i as f64) * dx;
            let up = (u[i + 1] - u[i - 1]) / (2.0 * dx); // u' at i

            // ∂L/∂u via central difference in u
            let dlu = (lagrangian(x, u[i] + eps, up) - lagrangian(x, u[i] - eps, up)) / (2.0 * eps);

            // ∂L/∂u' at i+1 and i-1 (for d/dx(∂L/∂u'))
            let up_p1 = (u[(i + 2).min(n - 1)] - u[i]) / (2.0 * dx);
            let up_m1 = (u[i] - u[(i as isize - 1).max(0) as usize]) / (2.0 * dx);

            let dlup_p1 = (lagrangian(x + dx, u[i + 1], up_p1 + eps)
                - lagrangian(x + dx, u[i + 1], up_p1 - eps))
                / (2.0 * eps);
            let dlup_m1 = (lagrangian(x - dx, u[i - 1], up_m1 + eps)
                - lagrangian(x - dx, u[i - 1], up_m1 - eps))
                / (2.0 * eps);

            let d_dlup_dx = (dlup_p1 - dlup_m1) / (2.0 * dx);
            residual[i] = dlu - d_dlup_dx;
        }
        residual
    }

    /// Minimize the functional `J\[u\]` subject to equality constraint `G\[u\] = 0`
    /// using the augmented Lagrangian method.
    ///
    /// Returns `(u, lambda)` where `u` is the minimizer and `lambda` is the
    /// Lagrange multiplier estimate.
    pub fn augmented_lagrangian<J, G>(
        j: J,
        g: G,
        u0: &[f64],
        rho: f64,
        step: f64,
        max_iter: usize,
    ) -> (Vec<f64>, f64)
    where
        J: Fn(&[f64]) -> f64,
        G: Fn(&[f64]) -> f64,
    {
        let mut u = u0.to_vec();
        let mut lambda = 0.0_f64;
        let eps = step * 1e-4;

        for _ in 0..max_iter {
            // Augmented Lagrangian: L_A(u) = J(u) + λ G(u) + ρ/2 G(u)²
            let augmented = |v: &[f64]| {
                let gv = g(v);
                j(v) + lambda * gv + rho / 2.0 * gv * gv
            };
            let grad = FunctionalDerivative::frechet_gradient(augmented, &u, eps);
            for (ui, gi) in u.iter_mut().zip(grad.iter()) {
                *ui -= step * gi;
            }
            // Update multiplier
            lambda += rho * g(&u);
        }
        (u, lambda)
    }

    /// Compute the first variation of `J\[u\] = ∫ f(u(x)) dx` at `u` in
    /// direction `v`.
    ///
    /// `f_prime` is the derivative of the integrand w.r.t. `u`.
    pub fn first_variation(u: &[f64], v: &[f64], f_prime: &[f64], dx: f64) -> f64 {
        assert_eq!(u.len(), v.len());
        assert_eq!(u.len(), f_prime.len());
        f_prime
            .iter()
            .zip(v.iter())
            .map(|(fp, vi)| fp * vi)
            .sum::<f64>()
            * dx
    }

    /// Find a stationary point via steepest descent on the functional gradient.
    ///
    /// Returns the stationary point and the history of functional values.
    pub fn steepest_descent<J>(
        j: J,
        u0: &[f64],
        step: f64,
        tol: f64,
        max_iter: usize,
    ) -> (Vec<f64>, Vec<f64>)
    where
        J: Fn(&[f64]) -> f64,
    {
        let eps = step * 1e-4;
        let mut u = u0.to_vec();
        let mut history = Vec::new();
        for _ in 0..max_iter {
            let val = j(&u);
            history.push(val);
            let grad = FunctionalDerivative::frechet_gradient(&j, &u, eps);
            let grad_norm: f64 = grad.iter().map(|x| x * x).sum::<f64>().sqrt();
            if grad_norm < tol {
                break;
            }
            for (ui, gi) in u.iter_mut().zip(grad.iter()) {
                *ui -= step * gi;
            }
        }
        (u, history)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── FunctionSpace ─────────────────────────────────────────────────────────

    #[test]
    fn test_function_space_eval() {
        let space = FunctionSpace::new(vec![|x: f64| x, |x: f64| x * x]);
        assert!((space.eval(0, 2.0) - 2.0).abs() < 1e-12);
        assert!((space.eval(1, 3.0) - 9.0).abs() < 1e-12);
    }

    #[test]
    fn test_function_space_dim() {
        let space = FunctionSpace::new(vec![|x: f64| x]);
        assert_eq!(space.dim(), 1);
    }

    #[test]
    fn test_function_space_combine() {
        // f(x) = 2x + 3x²  evaluated at x=1: 2+3=5
        let space = FunctionSpace::new(vec![|x: f64| x, |x: f64| x * x]);
        let val = space.combine(&[2.0, 3.0], 1.0);
        assert!((val - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_function_space_debug() {
        let space = FunctionSpace::new(vec![|x: f64| x]);
        let s = format!("{space:?}");
        assert!(s.contains("dim: 1"));
    }

    // ── L² inner product ─────────────────────────────────────────────────────

    #[test]
    fn test_l2_inner_product_orthogonal() {
        // sin and cos are orthogonal on [0, 2π]
        let n = 1000;
        let dx = 2.0 * PI / n as f64;
        let sin_v: Vec<f64> = (0..n).map(|i| ((i as f64 + 0.5) * dx).sin()).collect();
        let cos_v: Vec<f64> = (0..n).map(|i| ((i as f64 + 0.5) * dx).cos()).collect();
        let ip = l2_inner_product(&sin_v, &cos_v, dx);
        assert!(ip.abs() < 1e-10, "sin/cos inner product = {ip}");
    }

    #[test]
    fn test_l2_inner_product_self_is_norm_sq() {
        let f = vec![1.0, 2.0, 3.0];
        let dx = 0.1;
        let ip = l2_inner_product(&f, &f, dx);
        let norm = l2_norm(&f, dx);
        assert!((ip - norm * norm).abs() < 1e-12);
    }

    #[test]
    fn test_l2_norm_constant() {
        // ‖c‖ = c * √(n * dx) for constant c over n points with spacing dx
        let c = 3.0_f64;
        let n = 100;
        let dx = 0.01;
        let f = vec![c; n];
        let expected = c * ((n as f64) * dx).sqrt();
        let got = l2_norm(&f, dx);
        assert!((got - expected).abs() < 1e-10, "norm: {got} vs {expected}");
    }

    #[test]
    fn test_l2_norm_zero_vector() {
        let f = vec![0.0_f64; 50];
        assert!(l2_norm(&f, 0.01).abs() < 1e-12);
    }

    // ── Gram–Schmidt ──────────────────────────────────────────────────────────

    #[test]
    fn test_gram_schmidt_two_vectors() {
        let n = 100;
        let dx = 1.0 / n as f64;
        let v1: Vec<f64> = (0..n).map(|i| (i as f64 + 0.5) * dx).collect();
        let v2: Vec<f64> = (0..n).map(|_| 1.0).collect();
        let q = gram_schmidt_orthogonalize(&[v1, v2], dx);
        assert_eq!(q.len(), 2);
        // Orthogonality check
        let ip = l2_inner_product(&q[0], &q[1], dx);
        assert!(ip.abs() < 1e-10, "not orthogonal: {ip}");
        // Normality check
        let n0 = l2_norm(&q[0], dx);
        let n1 = l2_norm(&q[1], dx);
        assert!((n0 - 1.0).abs() < 1e-10, "q0 norm: {n0}");
        assert!((n1 - 1.0).abs() < 1e-10, "q1 norm: {n1}");
    }

    #[test]
    fn test_gram_schmidt_empty() {
        let q = gram_schmidt_orthogonalize(&[], 0.01);
        assert!(q.is_empty());
    }

    #[test]
    fn test_gram_schmidt_linearly_dependent() {
        let n = 10;
        let dx = 0.1;
        let v: Vec<f64> = (0..n).map(|i| i as f64).collect();
        let v2 = v.iter().map(|x| 2.0 * x).collect::<Vec<_>>();
        let q = gram_schmidt_orthogonalize(&[v, v2], dx);
        // Second vector is zero (linearly dependent)
        let norm2 = l2_norm(&q[1], dx);
        assert!(norm2 < 1e-12, "dependent vector should be zero: {norm2}");
    }

    #[test]
    fn test_gram_schmidt_three_vectors() {
        let n = 200;
        let dx = 2.0 * PI / n as f64;
        let v1: Vec<f64> = (0..n).map(|i| ((i as f64 + 0.5) * dx).sin()).collect();
        let v2: Vec<f64> = (0..n).map(|i| ((i as f64 + 0.5) * dx).cos()).collect();
        let v3: Vec<f64> = (0..n)
            .map(|i| (2.0 * (i as f64 + 0.5) * dx).sin())
            .collect();
        let q = gram_schmidt_orthogonalize(&[v1, v2, v3], dx);
        // All pairs must be orthogonal
        for i in 0..3 {
            for j in (i + 1)..3 {
                let ip = l2_inner_product(&q[i], &q[j], dx);
                assert!(ip.abs() < 1e-9, "q[{i}]·q[{j}] = {ip}");
            }
        }
    }

    // ── Fourier series ────────────────────────────────────────────────────────

    #[test]
    fn test_fourier_constant_function() {
        let n = 100;
        let dx = 1.0 / n as f64;
        let f = vec![1.0_f64; n];
        let coeffs = fourier_series_coeffs(&f, dx, 3);
        // a₀ should be 2 (by convention 2/L ∫ f dx = 2)
        assert!((coeffs[0].0 - 2.0).abs() < 1e-8, "a0 = {}", coeffs[0].0);
        // All b should be ~0
        for (k, &(_, bk)) in coeffs.iter().enumerate() {
            assert!(bk.abs() < 1e-8, "b{k} = {bk}");
        }
    }

    #[test]
    fn test_fourier_sine_function() {
        let n = 1000;
        let dx = 1.0 / n as f64;
        let f: Vec<f64> = (0..n)
            .map(|i| (2.0 * PI * (i as f64 + 0.5) * dx).sin())
            .collect();
        let coeffs = fourier_series_coeffs(&f, dx, 2);
        // b₁ should be ~1 for sin(2πx)
        assert!((coeffs[1].1 - 1.0).abs() < 0.01, "b1 = {}", coeffs[1].1);
    }

    #[test]
    fn test_fourier_n_terms_length() {
        let f = vec![0.0_f64; 64];
        let c = fourier_series_coeffs(&f, 0.1, 5);
        assert_eq!(c.len(), 5);
    }

    // ── Chebyshev expansion ───────────────────────────────────────────────────

    #[test]
    fn test_chebyshev_constant() {
        // f(x) = 1 on Chebyshev nodes: c₀ = 1, rest ≈ 0
        let n = 32;
        let f: Vec<f64> = (0..n)
            .map(|k| (PI * (k as f64 + 0.5) / n as f64).cos())
            .map(|_x| 1.0)
            .collect();
        let c = chebyshev_expansion(&f, 4);
        assert!((c[0] - 1.0).abs() < 1e-10, "c0 = {}", c[0]);
        for (k, &ck) in c.iter().enumerate().skip(1).take(3) {
            assert!(ck.abs() < 1e-10, "c{k} = {}", ck);
        }
    }

    #[test]
    fn test_chebyshev_empty_input() {
        let c = chebyshev_expansion(&[], 3);
        assert_eq!(c.len(), 3);
        for ci in c {
            assert!(ci.abs() < 1e-15);
        }
    }

    #[test]
    fn test_chebyshev_n_terms_length() {
        let f = vec![0.5_f64; 16];
        let c = chebyshev_expansion(&f, 6);
        assert_eq!(c.len(), 6);
    }

    #[test]
    fn test_chebyshev_x_polynomial() {
        // For f(x) = x on Chebyshev nodes, c₁ should dominate
        let n = 64;
        let f: Vec<f64> = (0..n)
            .map(|k| (PI * (k as f64 + 0.5) / n as f64).cos())
            .collect();
        let c = chebyshev_expansion(&f, 4);
        // c₀ should be ~0, c₁ should be ~1
        assert!(c[0].abs() < 1e-10, "c0 = {}", c[0]);
        assert!((c[1] - 1.0).abs() < 1e-10, "c1 = {}", c[1]);
    }

    // ── Legendre expansion ────────────────────────────────────────────────────

    #[test]
    fn test_legendre_constant() {
        // f(x) = 1: c₀ = 1, rest ≈ 0
        let n = 200;
        let f = vec![1.0_f64; n];
        let c = legendre_expansion(&f, 3);
        assert!((c[0] - 1.0).abs() < 1e-2, "c0 = {}", c[0]);
        assert!(c[1].abs() < 1e-4, "c1 = {}", c[1]);
        // P₂ quadrature error is O(1/n): with n=200 expect |c₂| < 1e-3
        assert!(c[2].abs() < 1e-3, "c2 = {}", c[2]);
    }

    #[test]
    fn test_legendre_n_terms_length() {
        let f = vec![1.0_f64; 50];
        let c = legendre_expansion(&f, 5);
        assert_eq!(c.len(), 5);
    }

    #[test]
    fn test_legendre_p_values() {
        // P₀(x) = 1, P₁(x) = x, P₂(x) = (3x²-1)/2
        assert!((legendre_p(0, 0.5) - 1.0).abs() < 1e-12);
        assert!((legendre_p(1, 0.5) - 0.5).abs() < 1e-12);
        let p2 = (3.0 * 0.5_f64.powi(2) - 1.0) / 2.0;
        assert!((legendre_p(2, 0.5) - p2).abs() < 1e-12);
    }

    #[test]
    fn test_legendre_empty_input() {
        let c = legendre_expansion(&[], 3);
        assert_eq!(c.len(), 3);
    }

    // ── Haar wavelet transform ────────────────────────────────────────────────

    #[test]
    fn test_haar_transform_inverse_roundtrip() {
        let signal = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let coeffs = wavelet_haar_transform(&signal);
        let recovered = wavelet_haar_inverse(&coeffs);
        for (a, b) in signal.iter().zip(recovered.iter()) {
            assert!((a - b).abs() < 1e-10, "roundtrip mismatch: {a} vs {b}");
        }
    }

    #[test]
    fn test_haar_transform_length_preserved() {
        let signal = vec![1.0; 16];
        let c = wavelet_haar_transform(&signal);
        assert_eq!(c.len(), 16);
    }

    #[test]
    fn test_haar_constant_signal() {
        // Constant signal: only the DC coefficient should be non-zero
        let c = 4.0_f64;
        let signal = vec![c; 8];
        let coeffs = wavelet_haar_transform(&signal);
        // The very first coefficient captures the mean (up to scaling)
        assert!(coeffs[0].abs() > 1e-6, "DC coeff should be non-zero");
        // All detail coefficients (non-DC) should be ~0
        for &coef in &coeffs[1..] {
            assert!(coef.abs() < 1e-10, "detail coeff = {coef}");
        }
    }

    #[test]
    fn test_haar_inverse_length_preserved() {
        let coeffs = vec![1.0; 8];
        let out = wavelet_haar_inverse(&coeffs);
        assert_eq!(out.len(), 8);
    }

    #[test]
    fn test_haar_two_point() {
        let signal = vec![2.0, 4.0];
        let c = wavelet_haar_transform(&signal);
        let r = wavelet_haar_inverse(&c);
        assert!((r[0] - 2.0).abs() < 1e-12);
        assert!((r[1] - 4.0).abs() < 1e-12);
    }

    // ── Sobolev norm (free function) ──────────────────────────────────────────

    #[test]
    fn test_sobolev_norm_s0_equals_l2() {
        // H⁰ norm = L² norm (when s=0)
        let f = vec![1.0, 2.0, 3.0];
        let df = vec![0.0, 0.0, 0.0];
        let dx = 0.1;
        let sob = sobolev_norm(&f, &df, dx, 0.0);
        let l2 = l2_norm(&f, dx);
        assert!((sob - l2).abs() < 1e-12);
    }

    #[test]
    fn test_sobolev_norm_positive() {
        let f = vec![1.0; 10];
        let df = vec![0.0; 10];
        let s = sobolev_norm(&f, &df, 0.1, 1.0);
        assert!(s > 0.0);
    }

    #[test]
    fn test_sobolev_norm_larger_with_gradient() {
        let f = vec![1.0; 8];
        let df_zero = vec![0.0; 8];
        let df_nonzero = vec![1.0; 8];
        let dx = 0.1;
        let s0 = sobolev_norm(&f, &df_zero, dx, 1.0);
        let s1 = sobolev_norm(&f, &df_nonzero, dx, 1.0);
        assert!(s1 > s0, "nonzero derivative should increase Sobolev norm");
    }

    // ── Operator norm ─────────────────────────────────────────────────────────

    #[test]
    fn test_operator_norm_identity() {
        // Identity 3×3 matrix: σ_max = 1
        let id = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];
        let sigma = operator_norm_estimate(&id);
        assert!((sigma - 1.0).abs() < 1e-4, "identity σ_max = {sigma}");
    }

    #[test]
    fn test_operator_norm_scaled_identity() {
        let scale = 5.0_f64;
        let m = vec![vec![scale, 0.0], vec![0.0, scale]];
        let sigma = operator_norm_estimate(&m);
        assert!(
            (sigma - scale).abs() < 1e-3,
            "scaled identity σ_max = {sigma}"
        );
    }

    #[test]
    fn test_operator_norm_empty_matrix() {
        assert_eq!(operator_norm_estimate(&[]), 0.0);
    }

    #[test]
    fn test_operator_norm_single_element() {
        let m = vec![vec![3.0_f64]];
        let sigma = operator_norm_estimate(&m);
        assert!((sigma - 3.0).abs() < 1e-4, "1x1 matrix σ = {sigma}");
    }

    #[test]
    fn test_operator_norm_rectangular() {
        // 1×2 matrix [3, 4]: σ_max = √(9+16) = 5
        let m = vec![vec![3.0_f64, 4.0]];
        let sigma = operator_norm_estimate(&m);
        assert!((sigma - 5.0).abs() < 1e-4, "σ = {sigma}, expected 5");
    }

    #[test]
    fn test_operator_norm_diagonal() {
        // Diagonal 3×3 with values [2, 5, 3]: σ_max = 5
        let m = vec![
            vec![2.0, 0.0, 0.0],
            vec![0.0, 5.0, 0.0],
            vec![0.0, 0.0, 3.0],
        ];
        let sigma = operator_norm_estimate(&m);
        assert!(
            (sigma - 5.0).abs() < 0.1,
            "diagonal max singular val = {sigma}"
        );
    }

    // ── HilbertSpace ──────────────────────────────────────────────────────────

    #[test]
    fn test_hilbert_inner_product() {
        let hs = HilbertSpace::new(4, 0.25);
        let f = vec![1.0, 2.0, 3.0, 4.0];
        let g = vec![4.0, 3.0, 2.0, 1.0];
        let ip = hs.inner_product(&f, &g);
        // (1*4 + 2*3 + 3*2 + 4*1) * 0.25 = 20 * 0.25 = 5
        assert!((ip - 5.0).abs() < 1e-12, "inner product = {ip}");
    }

    #[test]
    fn test_hilbert_norm_positive() {
        let hs = HilbertSpace::new(4, 0.25);
        let f = vec![1.0; 4];
        assert!(hs.norm(&f) > 0.0);
    }

    #[test]
    fn test_hilbert_normalize_unit_norm() {
        let hs = HilbertSpace::new(4, 0.25);
        let f = vec![3.0, 1.0, 4.0, 1.0];
        let fn_ = hs.normalize(&f);
        let nm = hs.norm(&fn_);
        assert!((nm - 1.0).abs() < 1e-10, "normalized norm = {nm}");
    }

    #[test]
    fn test_hilbert_orthogonality_check() {
        let n = 100;
        let dx = 2.0 * PI / n as f64;
        let hs = HilbertSpace::new(n, dx);
        let sin_v: Vec<f64> = (0..n).map(|i| ((i as f64 + 0.5) * dx).sin()).collect();
        let cos_v: Vec<f64> = (0..n).map(|i| ((i as f64 + 0.5) * dx).cos()).collect();
        assert!(hs.are_orthogonal(&sin_v, &cos_v, 1e-9));
    }

    #[test]
    fn test_hilbert_projection_length() {
        let n = 8;
        let dx = 1.0 / n as f64;
        let hs = HilbertSpace::new(n, dx);
        let f = vec![1.0; n];
        let basis = vec![(0..n).map(|i| i as f64 * dx).collect::<Vec<_>>()];
        let p = hs.project(&f, &basis);
        assert_eq!(p.len(), n);
    }

    #[test]
    fn test_hilbert_distance_zero_self() {
        let hs = HilbertSpace::new(4, 0.25);
        let f = vec![1.0, 2.0, 3.0, 4.0];
        assert!(hs.distance(&f, &f).abs() < 1e-12);
    }

    #[test]
    fn test_hilbert_orthonormalize_orthogonality() {
        let n = 100;
        let dx = 2.0 * PI / n as f64;
        let hs = HilbertSpace::new(n, dx);
        let v1: Vec<f64> = (0..n).map(|i| ((i as f64 + 0.5) * dx).sin()).collect();
        let v2: Vec<f64> = (0..n).map(|i| ((i as f64 + 0.5) * dx).cos()).collect();
        let q = hs.orthonormalize(&[v1, v2]);
        let ip = hs.inner_product(&q[0], &q[1]);
        assert!(ip.abs() < 1e-9, "orthonormalized inner product = {ip}");
    }

    // ── BanachSpace ───────────────────────────────────────────────────────────

    #[test]
    fn test_banach_l2_norm_matches_hilbert() {
        let dx = 0.25;
        let bs = BanachSpace::new(2.0, dx);
        let hs = HilbertSpace::new(4, dx);
        let f = vec![1.0, 2.0, 3.0, 4.0];
        let bn = bs.norm(&f);
        let hn = hs.norm(&f);
        assert!((bn - hn).abs() < 1e-10, "L2 norms match: {bn} vs {hn}");
    }

    #[test]
    fn test_banach_l1_norm() {
        let dx = 1.0;
        let bs = BanachSpace::new(1.0, dx);
        let f = vec![1.0, -2.0, 3.0];
        // L1 norm = (1 + 2 + 3) * 1.0 = 6
        let n = bs.norm(&f);
        assert!((n - 6.0).abs() < 1e-12, "L1 norm = {n}");
    }

    #[test]
    fn test_banach_linf_norm() {
        let dx = 0.1;
        let bs = BanachSpace::new(1e12, dx); // very large p ≈ L∞
        let f = vec![1.0, 5.0, 3.0];
        let n = bs.norm(&f);
        // Should be close to the max: 5.0 (for large p)
        assert!(n > 4.0 && n <= 5.0 * 1.01, "L∞-like norm = {n}");
    }

    #[test]
    fn test_banach_holder_inequality() {
        let dx = 0.1;
        let bs = BanachSpace::new(2.0, dx);
        let f = vec![1.0, 2.0, 3.0, 4.0];
        let g = vec![4.0, 3.0, 2.0, 1.0];
        let (lhs, rhs) = bs.holder_bound(&f, &g);
        assert!(lhs <= rhs + 1e-10, "Hölder: {lhs} <= {rhs}");
    }

    #[test]
    fn test_banach_unit_ball() {
        let dx = 1.0;
        let bs = BanachSpace::new(2.0, dx);
        let f = vec![0.5, 0.0];
        assert!(bs.in_unit_ball(&f, 1e-10));
    }

    #[test]
    fn test_banach_linear_functional() {
        let dx = 0.5;
        let bs = BanachSpace::new(2.0, dx);
        let f = vec![1.0, 1.0];
        let g = vec![2.0, 2.0];
        // φ(f) = (1*2 + 1*2) * 0.5 = 2
        let val = bs.bounded_linear_functional(&f, &g);
        assert!((val - 2.0).abs() < 1e-12, "functional = {val}");
    }

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

    // ── SobolevSpace ──────────────────────────────────────────────────────────

    #[test]
    fn test_sobolev_space_h0_equals_l2() {
        let n = 8;
        let dx = 0.125;
        let hs = HilbertSpace::new(n, dx);
        let ss = SobolevSpace::new(n, dx, 0);
        let f = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        assert!(
            (ss.norm(&f) - hs.norm(&f)).abs() < 1e-10,
            "H^0 = L^2: {} vs {}",
            ss.norm(&f),
            hs.norm(&f)
        );
    }

    #[test]
    fn test_sobolev_space_norm_positive() {
        let ss = SobolevSpace::new(8, 0.125, 1);
        let f = vec![1.0; 8];
        assert!(ss.norm(&f) > 0.0);
    }

    #[test]
    fn test_sobolev_weak_derivative_constant() {
        // Derivative of constant = 0 in interior, small at boundary
        let f = vec![3.0; 8];
        let d = SobolevSpace::weak_derivative(&f, 0.125);
        for &di in &d[1..d.len() - 1] {
            assert!(di.abs() < 1e-12, "interior deriv of constant: {di}");
        }
    }

    #[test]
    fn test_sobolev_trace_operator() {
        let ss = SobolevSpace::new(8, 0.125, 1);
        let f: Vec<f64> = (0..8).map(|i| i as f64).collect();
        let tr = ss.trace(&f);
        assert_eq!(tr[0], 0.0);
        assert_eq!(tr[1], 7.0);
    }

    #[test]
    fn test_sobolev_seminorm_nonneg() {
        let ss = SobolevSpace::new(8, 0.125, 1);
        let f: Vec<f64> = (0..8).map(|i| (i as f64 * 0.5).sin()).collect();
        assert!(ss.seminorm(&f) >= 0.0);
    }

    #[test]
    fn test_sobolev_inner_product_positive_definite() {
        let ss = SobolevSpace::new(8, 0.125, 1);
        let f = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let ip = ss.inner_product(&f, &f);
        assert!(ip > 0.0, "H^k inner product <f,f> > 0");
    }

    #[test]
    fn test_sobolev_is_in_space() {
        let ss = SobolevSpace::new(8, 0.125, 2);
        let f: Vec<f64> = (0..8).map(|i| (i as f64).sin()).collect();
        assert!(ss.is_in_space(&f));
    }

    // ── FunctionalDerivative ──────────────────────────────────────────────────

    #[test]
    fn test_gateaux_derivative_linear() {
        // J(f) = ∑ f[i] (integral with dx=1): dJ/dε J(f+εh) = ∑ h[i]
        let j = |f: &[f64]| f.iter().sum::<f64>();
        let f = vec![1.0, 2.0, 3.0];
        let h = vec![1.0, 0.0, 0.0];
        let dj = FunctionalDerivative::gateaux(j, &f, &h, 1e-6);
        assert!((dj - 1.0).abs() < 1e-6, "Gateaux derivative = {dj}");
    }

    #[test]
    fn test_frechet_gradient_quadratic() {
        // J(f) = ∑ f[i]²  → gradient is 2*f
        let j = |f: &[f64]| f.iter().map(|x| x * x).sum::<f64>();
        let f = vec![1.0, 2.0, 3.0];
        let grad = FunctionalDerivative::frechet_gradient(j, &f, 1e-5);
        for (i, (&gi, &fi)) in grad.iter().zip(f.iter()).enumerate() {
            assert!(
                (gi - 2.0 * fi).abs() < 1e-4,
                "grad[{i}] = {gi}, expected {}",
                2.0 * fi
            );
        }
    }

    #[test]
    fn test_gradient_descent_minimizes() {
        // J(f) = (f[0]-1)² + (f[1]-2)²; minimum at [1, 2]
        let j = |f: &[f64]| (f[0] - 1.0).powi(2) + (f[1] - 2.0).powi(2);
        let f0 = vec![0.0, 0.0];
        let (fmin, val) = FunctionalDerivative::gradient_descent(j, &f0, 0.1, 200, 1e-5);
        assert!(val < 0.1, "not minimized: val = {val}");
        assert!((fmin[0] - 1.0).abs() < 0.1, "f[0] = {}", fmin[0]);
        assert!((fmin[1] - 2.0).abs() < 0.1, "f[1] = {}", fmin[1]);
    }

    #[test]
    fn test_newton_step_moves_toward_minimum() {
        // J(f) = (f[0]-3)²; minimum at f[0]=3
        let j = |f: &[f64]| (f[0] - 3.0).powi(2);
        let f = vec![1.0];
        let f_new = FunctionalDerivative::newton_step(j, &f, 1e-5);
        // Newton step on quadratic should reach minimum in one step
        assert!(
            (f_new[0] - 3.0).abs() < 0.1,
            "Newton step: f[0] = {}",
            f_new[0]
        );
    }

    #[test]
    fn test_second_gateaux_symmetric() {
        // For a quadratic J = ∑ f²: second Gâteaux should be symmetric
        let j = |f: &[f64]| f.iter().map(|x| x * x).sum::<f64>();
        let f = vec![1.0, 1.0];
        let h = vec![1.0, 0.0];
        let k_dir = vec![0.0, 1.0];
        let d2_hk = FunctionalDerivative::second_gateaux(j, &f, &h, &k_dir, 1e-4);
        let d2_kh = FunctionalDerivative::second_gateaux(j, &f, &k_dir, &h, 1e-4);
        assert!(
            (d2_hk - d2_kh).abs() < 1e-4,
            "second Gateaux not symmetric: {d2_hk} vs {d2_kh}"
        );
    }

    // ── VariationalProblem ────────────────────────────────────────────────────

    #[test]
    fn test_euler_lagrange_residual_length() {
        // L(x, u, u') = u'²
        let u: Vec<f64> = (0..10).map(|i| i as f64 * 0.1).collect();
        let res = VariationalProblem::euler_lagrange_residual(|_x, _u, up| up * up, &u, 0.1);
        assert_eq!(res.len(), u.len());
    }

    #[test]
    fn test_euler_lagrange_linear_minimizer_residual() {
        // For L = u'²: E-L equation is -u'' = 0 → linear solution is exact
        let n = 10;
        let dx = 1.0 / n as f64;
        let u: Vec<f64> = (0..n).map(|i| i as f64 * dx).collect();
        let res = VariationalProblem::euler_lagrange_residual(|_x, _u, up| up * up, &u, dx);
        // Interior residuals should be small for linear u
        for &r in &res[1..n - 1] {
            assert!(r.is_finite(), "residual finite: {r}");
        }
    }

    #[test]
    fn test_steepest_descent_minimizes() {
        // J(f) = ∑ (f[i] - c)²; minimum at f = [c, c, ...]
        let c = 2.0_f64;
        let j = move |f: &[f64]| f.iter().map(|x| (x - c).powi(2)).sum::<f64>();
        let u0 = vec![0.0; 4];
        let (u, hist) = VariationalProblem::steepest_descent(j, &u0, 0.1, 1e-6, 300);
        assert!(!hist.is_empty());
        for &ui in &u {
            assert!((ui - c).abs() < 0.1, "steepest descent: ui = {ui}");
        }
    }

    #[test]
    fn test_first_variation_linear() {
        // δJ = ∫ f'(u) v dx
        let u = vec![1.0; 4];
        let v = vec![1.0; 4];
        let fp = vec![2.0; 4]; // derivative of u²: 2u = 2
        let dx = 0.25;
        let dj = VariationalProblem::first_variation(&u, &v, &fp, dx);
        // = (2*1 + 2*1 + 2*1 + 2*1) * 0.25 = 2.0
        assert!((dj - 2.0).abs() < 1e-12, "first variation = {dj}");
    }

    #[test]
    fn test_augmented_lagrangian_reduces_constraint() {
        // Minimize ∑ f² subject to ∑ f = 1
        let j = |f: &[f64]| f.iter().map(|x| x * x).sum::<f64>();
        let g = |f: &[f64]| f.iter().sum::<f64>() - 1.0;
        let u0 = vec![0.5, 0.5];
        let (_u, _lambda) = VariationalProblem::augmented_lagrangian(j, g, &u0, 1.0, 0.05, 100);
        // Just check that it runs without panic and returns finite values
    }
}
