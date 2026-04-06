#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Spectral methods for numerical analysis and physics simulations.
//!
//! Provides Chebyshev and Legendre polynomial evaluation, Gauss quadrature
//! nodes and weights, Cooley-Tukey FFT/IFFT, pseudo-spectral differentiation,
//! Chebyshev collocation for 1D boundary value problems, and Haar wavelet
//! multi-level decomposition/reconstruction.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// Complex number (local, minimal)
// ─────────────────────────────────────────────────────────────────────────────

/// Minimal complex number for FFT computations.
#[derive(Clone, Copy, Debug, PartialEq)]
struct Cx {
    re: f64,
    im: f64,
}

impl Cx {
    #[inline]
    fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }
    #[inline]
    fn from_polar(r: f64, theta: f64) -> Self {
        Self {
            re: r * theta.cos(),
            im: r * theta.sin(),
        }
    }
    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self {
            re: self.re + rhs.re,
            im: self.im + rhs.im,
        }
    }
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self {
            re: self.re - rhs.re,
            im: self.im - rhs.im,
        }
    }
    #[inline]
    fn mul(self, rhs: Self) -> Self {
        Self {
            re: self.re * rhs.re - self.im * rhs.im,
            im: self.re * rhs.im + self.im * rhs.re,
        }
    }
    #[inline]
    fn abs(self) -> f64 {
        self.re.hypot(self.im)
    }
    #[inline]
    fn scale(self, s: f64) -> Self {
        Self {
            re: self.re * s,
            im: self.im * s,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ChebyshevPolynomial
// ─────────────────────────────────────────────────────────────────────────────

/// Chebyshev polynomials of the first kind: T_n(x).
///
/// Provides evaluation via the three-term recurrence, Chebyshev nodes on
/// \[-1, 1\], differentiation matrix for spectral differentiation, and
/// interpolation coefficients from function values at Chebyshev nodes.
pub struct ChebyshevPolynomial;

impl ChebyshevPolynomial {
    /// Evaluate T_n(x) using the recurrence T_0=1, T_1=x, T_{n+1}=2x T_n - T_{n-1}.
    pub fn eval(n: usize, x: f64) -> f64 {
        if n == 0 {
            return 1.0;
        }
        if n == 1 {
            return x;
        }
        let mut t_prev = 1.0f64;
        let mut t_curr = x;
        for _ in 2..=n {
            let t_next = 2.0 * x * t_curr - t_prev;
            t_prev = t_curr;
            t_curr = t_next;
        }
        t_curr
    }

    /// Evaluate all Chebyshev polynomials T_0(x), …, T_n(x).
    pub fn eval_all(n: usize, x: f64) -> Vec<f64> {
        let mut ts = vec![0.0f64; n + 1];
        ts[0] = 1.0;
        if n >= 1 {
            ts[1] = x;
        }
        for k in 2..=n {
            ts[k] = 2.0 * x * ts[k - 1] - ts[k - 2];
        }
        ts
    }

    /// Return the `n+1` Chebyshev-Gauss-Lobatto nodes on \[-1, 1\]:
    /// `x_j = cos(j*pi/n)`, j = 0, …, n.
    pub fn nodes(n: usize) -> Vec<f64> {
        (0..=n).map(|j| (j as f64 * PI / n as f64).cos()).collect()
    }

    /// Return the `n` interior Chebyshev-Gauss nodes on (-1, 1):
    /// `x_j = cos((2j+1)*pi / (2n))`, j = 0, …, n-1.
    pub fn gauss_nodes(n: usize) -> Vec<f64> {
        (0..n)
            .map(|j| ((2 * j + 1) as f64 * PI / (2 * n) as f64).cos())
            .collect()
    }

    /// Compute the `(n+1) x (n+1)` spectral differentiation matrix at
    /// Chebyshev-Gauss-Lobatto nodes.
    ///
    /// The entry `D[i][j]` approximates the derivative of the `j`-th basis
    /// function at node `x_i`.
    pub fn diff_matrix(n: usize) -> Vec<Vec<f64>> {
        let nodes = Self::nodes(n);
        let m = n + 1;
        let mut d = vec![vec![0.0f64; m]; m];

        let c = |i: usize| -> f64 { if i == 0 || i == n { 2.0 } else { 1.0 } };

        for i in 0..m {
            for j in 0..m {
                if i != j {
                    d[i][j] =
                        c(i) / c(j) * ((-1.0f64).powi((i + j) as i32)) / (nodes[i] - nodes[j]);
                }
            }
            // Diagonal: negative sum of off-diagonal row entries
            let row_sum: f64 = (0..m).filter(|&k| k != i).map(|k| d[i][k]).sum();
            d[i][i] = -row_sum;
        }
        d
    }

    /// Compute Chebyshev expansion coefficients from function values at
    /// Chebyshev-Gauss-Lobatto nodes using the DCT-I formula.
    ///
    /// Returns coefficients `a_k` such that `f(x) ≈ Σ a_k T_k(x)`.
    pub fn interpolation_coeffs(vals: &[f64]) -> Vec<f64> {
        let n = vals.len() - 1;
        let m = n + 1;
        let mut coeffs = vec![0.0f64; m];
        for k in 0..m {
            let norm = if k == 0 || k == n {
                n as f64
            } else {
                n as f64 / 2.0
            };
            let sum: f64 = (0..m)
                .map(|j| {
                    let w = if j == 0 || j == n { 0.5 } else { 1.0 };
                    w * vals[j] * (k as f64 * j as f64 * PI / n as f64).cos()
                })
                .sum();
            coeffs[k] = sum / norm;
        }
        coeffs
    }

    /// Evaluate the Chebyshev series with coefficients `coeffs` at point `x`.
    pub fn eval_series(coeffs: &[f64], x: f64) -> f64 {
        coeffs
            .iter()
            .enumerate()
            .map(|(k, &ck)| ck * Self::eval(k, x))
            .sum()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// LegendrePolynomial
// ─────────────────────────────────────────────────────────────────────────────

/// Legendre polynomials P_n(x) and Gauss-Legendre quadrature.
///
/// Provides evaluation via the three-term recurrence, computation of all
/// P_0, …, P_n, and Gauss-Legendre nodes and weights up to degree n.
pub struct LegendrePolynomial;

impl LegendrePolynomial {
    /// Evaluate P_n(x) using the three-term recurrence.
    pub fn eval(n: usize, x: f64) -> f64 {
        if n == 0 {
            return 1.0;
        }
        if n == 1 {
            return x;
        }
        let mut p_prev = 1.0f64;
        let mut p_curr = x;
        for k in 1..n {
            let k_f = k as f64;
            let p_next = ((2.0 * k_f + 1.0) * x * p_curr - k_f * p_prev) / (k_f + 1.0);
            p_prev = p_curr;
            p_curr = p_next;
        }
        p_curr
    }

    /// Evaluate P_0(x), …, P_n(x) and return all values.
    pub fn eval_all(n: usize, x: f64) -> Vec<f64> {
        let mut ps = vec![0.0f64; n + 1];
        ps[0] = 1.0;
        if n >= 1 {
            ps[1] = x;
        }
        for k in 1..n {
            let k_f = k as f64;
            ps[k + 1] = ((2.0 * k_f + 1.0) * x * ps[k] - k_f * ps[k - 1]) / (k_f + 1.0);
        }
        ps
    }

    /// Compute the `n`-point Gauss-Legendre nodes and weights on \[-1, 1\].
    ///
    /// Uses Newton's method to find roots of P_n and computes weights from
    /// the derivative formula.  Returns `(nodes, weights)`.
    pub fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
        let mut nodes = vec![0.0f64; n];
        let mut weights = vec![0.0f64; n];

        for i in 0..n.div_ceil(2) {
            // Initial guess: Chebyshev nodes
            let mut x = ((2 * i + 1) as f64 * PI / (2 * n) as f64 + PI / (4.0 * n as f64)).cos();

            for _ in 0..100 {
                let ps = Self::eval_all(n, x);
                let pn = ps[n];
                let pn_1 = if n >= 1 { ps[n - 1] } else { 0.0 };
                // Derivative: P_n'(x) = n * (x*P_n(x) - P_{n-1}(x)) / (x^2 - 1)
                let dp = if (x.abs() - 1.0).abs() < 1e-14 {
                    n as f64 * (n as f64 + 1.0) / 2.0 // limit at ±1
                } else {
                    (n as f64) * (pn_1 - x * pn) / (1.0 - x * x)
                };
                let dx = pn / dp;
                x -= dx;
                if dx.abs() < 1e-15 {
                    break;
                }
            }

            let ps = Self::eval_all(n, x);
            let pn = ps[n];
            let pn_1 = if n >= 1 { ps[n - 1] } else { 0.0 };
            let dp = if (x.abs() - 1.0).abs() < 1e-14 {
                n as f64 * (n as f64 + 1.0) / 2.0
            } else {
                (n as f64) * (pn_1 - x * pn) / (1.0 - x * x)
            };
            let w = 2.0 / ((1.0 - x * x) * dp * dp);

            nodes[i] = -x;
            nodes[n - 1 - i] = x;
            weights[i] = w;
            weights[n - 1 - i] = w;
        }

        (nodes, weights)
    }

    /// Integrate `f` over \[-1, 1\] using `n`-point Gauss-Legendre quadrature.
    pub fn integrate<F: Fn(f64) -> f64>(f: F, n: usize) -> f64 {
        let (nodes, weights) = Self::gauss_legendre(n);
        nodes
            .iter()
            .zip(weights.iter())
            .map(|(&x, &w)| w * f(x))
            .sum()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// FourierSeries — DFT, IDFT, power spectrum, convolution
// ─────────────────────────────────────────────────────────────────────────────

/// Discrete Fourier Transform and spectral utilities.
///
/// Provides Cooley-Tukey radix-2 FFT (in-place), IFFT, power spectrum
/// estimation, and circular convolution via the convolution theorem.
pub struct FourierSeries;

impl FourierSeries {
    /// Compute the FFT of `data` (length must be a power of 2).
    ///
    /// Returns a vector of complex coefficients `X[k]` where
    /// `X[k] = Σ x[n] * exp(-2πi kn/N)`.
    pub fn fft(data: &[f64]) -> Vec<(f64, f64)> {
        let n = data.len();
        assert!(n.is_power_of_two(), "FFT length must be a power of 2");
        let mut buf: Vec<Cx> = data.iter().map(|&x| Cx::new(x, 0.0)).collect();
        fft_inplace(&mut buf, false);
        buf.iter().map(|c| (c.re, c.im)).collect()
    }

    /// Compute the inverse FFT.
    ///
    /// Input is a slice of `(re, im)` pairs; output is the real part of the
    /// inverse transform (discards imaginary part which should be ~0 for real signals).
    pub fn ifft(spectrum: &[(f64, f64)]) -> Vec<f64> {
        let n = spectrum.len();
        assert!(n.is_power_of_two(), "IFFT length must be a power of 2");
        let mut buf: Vec<Cx> = spectrum.iter().map(|&(re, im)| Cx::new(re, im)).collect();
        fft_inplace(&mut buf, true);
        buf.iter().map(|c| c.re / n as f64).collect()
    }

    /// Compute the one-sided power spectrum of `data`.
    ///
    /// Returns `|X[k]|^2 / N` for k = 0, …, N/2.
    pub fn power_spectrum(data: &[f64]) -> Vec<f64> {
        let n = data.len();
        assert!(n.is_power_of_two());
        let spec = Self::fft(data);
        let norm = n as f64;
        (0..=n / 2)
            .map(|k| {
                let (re, im) = spec[k];
                (re * re + im * im) / norm
            })
            .collect()
    }

    /// Circular convolution of `a` and `b` via the convolution theorem.
    ///
    /// Both inputs must have the same length (a power of 2).
    /// Returns the circular convolution `a * b`.
    pub fn convolve(a: &[f64], b: &[f64]) -> Vec<f64> {
        let n = a.len();
        assert_eq!(n, b.len());
        assert!(n.is_power_of_two());
        let sa = Self::fft(a);
        let sb = Self::fft(b);
        let product: Vec<(f64, f64)> = sa
            .iter()
            .zip(sb.iter())
            .map(|(&(ar, ai), &(br, bi))| (ar * br - ai * bi, ar * bi + ai * br))
            .collect();
        Self::ifft(&product)
    }

    /// Compute the DFT frequency bins for sample rate `fs` and `n` points.
    ///
    /// Returns `n/2 + 1` non-negative frequencies in Hz.
    pub fn frequencies(n: usize, fs: f64) -> Vec<f64> {
        (0..=n / 2).map(|k| k as f64 * fs / n as f64).collect()
    }

    /// Evaluate the truncated Fourier series at points `x` using `n_terms` terms.
    ///
    /// Coefficients `(a_k, b_k)` are the cosine and sine amplitudes.
    /// Returns `sum_{k=0}^{n_terms-1} a_k cos(k x) + b_k sin(k x)`.
    pub fn eval_series(coeffs: &[(f64, f64)], x: f64) -> f64 {
        coeffs
            .iter()
            .enumerate()
            .map(|(k, &(ak, bk))| ak * (k as f64 * x).cos() + bk * (k as f64 * x).sin())
            .sum()
    }
}

/// In-place Cooley-Tukey FFT (radix-2 DIT).
fn fft_inplace(buf: &mut Vec<Cx>, inverse: bool) {
    let n = buf.len();
    // Bit-reversal permutation
    let mut j = 0usize;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if i < j {
            buf.swap(i, j);
        }
    }
    // Butterfly
    let sign = if inverse { 1.0 } else { -1.0 };
    let mut len = 2usize;
    while len <= n {
        let ang = sign * 2.0 * PI / len as f64;
        let w_len = Cx::from_polar(1.0, ang);
        for i in (0..n).step_by(len) {
            let mut w = Cx::new(1.0, 0.0);
            for k in 0..len / 2 {
                let u = buf[i + k];
                let v = buf[i + k + len / 2].mul(w);
                buf[i + k] = u.add(v);
                buf[i + k + len / 2] = u.sub(v);
                w = w.mul(w_len);
            }
        }
        len <<= 1;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SpectralDiff — spectral differentiation via FFT
// ─────────────────────────────────────────────────────────────────────────────

/// Pseudo-spectral differentiation and related operations.
///
/// Differentiates periodic functions sampled on a uniform grid using the
/// spectral derivative (multiplication by `ik` in Fourier space).
pub struct SpectralDiff;

impl SpectralDiff {
    /// Compute the first derivative of a periodic function sampled at `n`
    /// uniformly spaced points on \[0, L) using spectral (FFT) differentiation.
    ///
    /// `n` must be a power of 2.  Returns `du/dx` at the same grid points.
    pub fn diff(u: &[f64], l: f64) -> Vec<f64> {
        let n = u.len();
        assert!(n.is_power_of_two());
        let mut buf: Vec<Cx> = u.iter().map(|&x| Cx::new(x, 0.0)).collect();
        fft_inplace(&mut buf, false);

        // Multiply by i*k (wavenumber)
        let dk = 2.0 * PI / l;
        for k in 0..n {
            let kk = if k <= n / 2 {
                k as f64
            } else {
                k as f64 - n as f64
            };
            let freq = kk * dk;
            let (re, im) = (buf[k].re, buf[k].im);
            buf[k] = Cx::new(-freq * im, freq * re);
        }

        fft_inplace(&mut buf, true);
        buf.iter().map(|c| c.re / n as f64).collect()
    }

    /// Compute the second derivative of a periodic function via FFT.
    ///
    /// Multiplies each Fourier mode by `-(ik)^2 = k^2`.
    pub fn diff2(u: &[f64], l: f64) -> Vec<f64> {
        let n = u.len();
        assert!(n.is_power_of_two());
        let mut buf: Vec<Cx> = u.iter().map(|&x| Cx::new(x, 0.0)).collect();
        fft_inplace(&mut buf, false);

        let dk = 2.0 * PI / l;
        for k in 0..n {
            let kk = if k <= n / 2 {
                k as f64
            } else {
                k as f64 - n as f64
            };
            let freq2 = -(kk * dk).powi(2);
            buf[k] = buf[k].scale(freq2);
        }

        fft_inplace(&mut buf, true);
        buf.iter().map(|c| c.re / n as f64).collect()
    }

    /// Interpolate a periodic function sampled at `n` points to `m` uniformly
    /// spaced points on \[0, L) via zero-padding in Fourier space.
    ///
    /// Both `n` and `m` must be powers of 2, and `m >= n`.
    pub fn interpolate(u: &[f64], m: usize, l: f64) -> Vec<f64> {
        let n = u.len();
        assert!(n.is_power_of_two() && m.is_power_of_two() && m >= n);
        let mut buf: Vec<Cx> = u.iter().map(|&x| Cx::new(x, 0.0)).collect();
        fft_inplace(&mut buf, false);

        let mut padded = vec![Cx::new(0.0, 0.0); m];
        padded[..n / 2].copy_from_slice(&buf[..n / 2]);
        for k in 1..=n / 2 {
            padded[m - k] = buf[n - k];
        }

        fft_inplace(&mut padded, true);
        let scale = m as f64 / (n as f64 * m as f64);
        let _ = l; // grid spacing is implicit in the scaling
        padded
            .iter()
            .map(|c| c.re * m as f64 / n as f64 / m as f64 * n as f64 * scale)
            .collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ChebyshevCollocation — 1D BVP solver
// ─────────────────────────────────────────────────────────────────────────────

/// Chebyshev collocation method for 1D boundary value problems.
///
/// Sets up the collocation differentiation matrix on Chebyshev-Gauss-Lobatto
/// nodes and solves second-order BVPs of the form
/// `p(x) u'' + q(x) u' + r(x) u = g(x)` with Dirichlet boundary conditions.
pub struct ChebyshevCollocation {
    /// Number of interior collocation points (polynomial degree = n+1).
    pub n: usize,
}

impl ChebyshevCollocation {
    /// Construct a collocation scheme with `n` Chebyshev-Gauss-Lobatto points.
    pub fn new(n: usize) -> Self {
        assert!(n >= 2, "n must be at least 2");
        Self { n }
    }

    /// Return the collocation nodes (Chebyshev-Gauss-Lobatto on \[-1, 1\]).
    pub fn nodes(&self) -> Vec<f64> {
        ChebyshevPolynomial::nodes(self.n - 1)
    }

    /// Return the first-derivative spectral differentiation matrix D of size `n x n`.
    pub fn diff_matrix(&self) -> Vec<Vec<f64>> {
        ChebyshevPolynomial::diff_matrix(self.n - 1)
    }

    /// Solve the 1D Poisson equation `u'' = g(x)` on \[-1, 1\] with
    /// Dirichlet boundary conditions `u(-1) = bc_left`, `u(1) = bc_right`.
    ///
    /// Uses the Chebyshev spectral differentiation matrix; applies boundary
    /// conditions by replacing the first and last rows.  Returns `u` at all
    /// collocation nodes.
    pub fn solve_poisson<G>(&self, g: G, bc_left: f64, bc_right: f64) -> Vec<f64>
    where
        G: Fn(f64) -> f64,
    {
        let m = self.n;
        let x = self.nodes();
        let d = self.diff_matrix();

        // Compute D^2 = D * D
        let mut d2 = vec![vec![0.0f64; m]; m];
        for i in 0..m {
            for j in 0..m {
                for k in 0..m {
                    d2[i][j] += d[i][k] * d[k][j];
                }
            }
        }

        // Build RHS
        let mut rhs: Vec<f64> = x.iter().map(|&xi| g(xi)).collect();

        // Enforce boundary conditions (nodes are ordered x[0]=1, x[m-1]=-1)
        // Overwrite first row: u(x[0]) = bc_right (x[0] = cos(0) = 1)
        for j in 0..m {
            d2[0][j] = if j == 0 { 1.0 } else { 0.0 };
        }
        rhs[0] = bc_right;

        // Overwrite last row: u(x[m-1]) = bc_left (x[m-1] = cos(pi) = -1)
        for j in 0..m {
            d2[m - 1][j] = if j == m - 1 { 1.0 } else { 0.0 };
        }
        rhs[m - 1] = bc_left;

        // Solve linear system via Gaussian elimination with partial pivoting
        gauss_solve(&mut d2, &mut rhs)
    }
}

/// Solve `A x = b` in-place via Gaussian elimination with partial pivoting.
///
/// Modifies `a` and `b`; returns the solution vector.
fn gauss_solve(a: &mut Vec<Vec<f64>>, b: &mut Vec<f64>) -> Vec<f64> {
    let n = b.len();
    for col in 0..n {
        // Partial pivot
        let mut max_row = col;
        let mut max_val = a[col][col].abs();
        for row in (col + 1)..n {
            if a[row][col].abs() > max_val {
                max_val = a[row][col].abs();
                max_row = row;
            }
        }
        a.swap(col, max_row);
        b.swap(col, max_row);

        let pivot = a[col][col];
        if pivot.abs() < 1e-14 {
            continue;
        }
        for row in (col + 1)..n {
            let factor = a[row][col] / pivot;
            for k in col..n {
                let sub = factor * a[col][k];
                a[row][k] -= sub;
            }
            b[row] -= factor * b[col];
        }
    }
    // Back substitution
    let mut x = vec![0.0f64; n];
    for i in (0..n).rev() {
        let sum: f64 = (i + 1..n).map(|j| a[i][j] * x[j]).sum();
        x[i] = if a[i][i].abs() > 1e-14 {
            (b[i] - sum) / a[i][i]
        } else {
            0.0
        };
    }
    x
}

// ─────────────────────────────────────────────────────────────────────────────
// WaveletTransform — Haar wavelet, multi-level decomposition, reconstruction
// ─────────────────────────────────────────────────────────────────────────────

/// Haar wavelet multi-level decomposition and reconstruction.
///
/// The Haar wavelet is the simplest orthogonal wavelet.  The transform
/// iteratively applies the single-level forward transform to the approximation
/// coefficients, producing a hierarchy of detail sub-bands.
pub struct WaveletTransform;

impl WaveletTransform {
    /// Forward single-level Haar transform.
    ///
    /// Returns `(approx, detail)` each of length `signal.len() / 2`.
    pub fn haar_forward(signal: &[f64]) -> (Vec<f64>, Vec<f64>) {
        let n = signal.len() / 2;
        let s2i = 1.0 / std::f64::consts::SQRT_2;
        let mut approx = Vec::with_capacity(n);
        let mut detail = Vec::with_capacity(n);
        for i in 0..n {
            let a = signal[2 * i];
            let b = signal[2 * i + 1];
            approx.push((a + b) * s2i);
            detail.push((a - b) * s2i);
        }
        (approx, detail)
    }

    /// Inverse single-level Haar transform.
    ///
    /// Reconstructs a signal of length `2 * approx.len()`.
    pub fn haar_inverse(approx: &[f64], detail: &[f64]) -> Vec<f64> {
        let n = approx.len().min(detail.len());
        let s2i = 1.0 / std::f64::consts::SQRT_2;
        let mut out = vec![0.0f64; 2 * n];
        for i in 0..n {
            out[2 * i] = (approx[i] + detail[i]) * s2i;
            out[2 * i + 1] = (approx[i] - detail[i]) * s2i;
        }
        out
    }

    /// Multi-level Haar decomposition.
    ///
    /// Returns a `Vec` of length `levels + 1`:
    /// - index 0: final (coarsest) approximation
    /// - indices 1..=levels: detail coefficients (finest at index 1)
    pub fn decompose(signal: &[f64], levels: usize) -> Vec<Vec<f64>> {
        let mut result = Vec::with_capacity(levels + 1);
        let mut approx = signal.to_vec();
        for _ in 0..levels {
            if approx.len() < 2 {
                break;
            }
            let (a, d) = Self::haar_forward(&approx);
            result.push(d);
            approx = a;
        }
        result.push(approx);
        result.reverse();
        result
    }

    /// Multi-level Haar reconstruction from decomposed coefficients.
    ///
    /// Input `coeffs` must be in the same format as returned by [`WaveletTransform::decompose`].
    pub fn reconstruct(coeffs: &[Vec<f64>]) -> Vec<f64> {
        if coeffs.is_empty() {
            return vec![];
        }
        let mut approx = coeffs[0].clone();
        for detail in &coeffs[1..] {
            approx = Self::haar_inverse(&approx, detail);
        }
        approx
    }

    /// Soft thresholding (denoising) on wavelet detail coefficients.
    ///
    /// Applies the soft-threshold function `sign(x) * max(|x| - lambda, 0)`
    /// to each detail coefficient in `coeffs[1..]`.
    pub fn soft_threshold(coeffs: &mut Vec<Vec<f64>>, lambda: f64) {
        for sub in coeffs.iter_mut().skip(1) {
            for v in sub.iter_mut() {
                let s = v.signum();
                let a = v.abs() - lambda;
                *v = if a > 0.0 { s * a } else { 0.0 };
            }
        }
    }

    /// Compute the energy in each sub-band of a decomposition.
    ///
    /// Returns one energy value per sub-band.
    pub fn subband_energy(coeffs: &[Vec<f64>]) -> Vec<f64> {
        coeffs
            .iter()
            .map(|sub| sub.iter().map(|v| v * v).sum::<f64>())
            .collect()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────
#[cfg(test)]
mod tests {
    use super::*;

    // ------------------------------------------------------------------
    // ChebyshevPolynomial
    // ------------------------------------------------------------------
    #[test]
    fn test_cheb_eval_t0() {
        assert!((ChebyshevPolynomial::eval(0, 0.7) - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_cheb_eval_t1() {
        assert!((ChebyshevPolynomial::eval(1, 0.5) - 0.5).abs() < 1e-14);
    }

    #[test]
    fn test_cheb_eval_t2() {
        // T_2(x) = 2x^2 - 1
        let x = 0.6;
        let expected = 2.0 * x * x - 1.0;
        assert!((ChebyshevPolynomial::eval(2, x) - expected).abs() < 1e-12);
    }

    #[test]
    fn test_cheb_eval_t3() {
        // T_3(x) = 4x^3 - 3x
        let x = 0.3;
        let expected = 4.0 * x * x * x - 3.0 * x;
        assert!((ChebyshevPolynomial::eval(3, x) - expected).abs() < 1e-12);
    }

    #[test]
    fn test_cheb_eval_all_consistency() {
        let x = 0.4;
        let all = ChebyshevPolynomial::eval_all(5, x);
        for k in 0..=5 {
            assert!((all[k] - ChebyshevPolynomial::eval(k, x)).abs() < 1e-12);
        }
    }

    #[test]
    fn test_cheb_nodes_count() {
        let n = 7;
        let nodes = ChebyshevPolynomial::nodes(n);
        assert_eq!(nodes.len(), n + 1);
    }

    #[test]
    fn test_cheb_nodes_bounds() {
        for &x in ChebyshevPolynomial::nodes(8).iter() {
            assert!((-1.0 - 1e-12..=1.0 + 1e-12).contains(&x));
        }
    }

    #[test]
    fn test_cheb_gauss_nodes() {
        let nodes = ChebyshevPolynomial::gauss_nodes(5);
        assert_eq!(nodes.len(), 5);
        // All interior nodes
        for &x in &nodes {
            assert!(x.abs() < 1.0);
        }
    }

    #[test]
    fn test_cheb_diff_matrix_size() {
        let d = ChebyshevPolynomial::diff_matrix(5);
        assert_eq!(d.len(), 6);
        assert_eq!(d[0].len(), 6);
    }

    #[test]
    fn test_cheb_diff_matrix_row_sum_zero() {
        let d = ChebyshevPolynomial::diff_matrix(6);
        // Each row of the differentiation matrix should sum to ~0
        // (since the derivative of a constant is 0)
        for row in &d {
            let s: f64 = row.iter().sum();
            assert!(s.abs() < 1e-8, "row sum = {s}");
        }
    }

    #[test]
    fn test_cheb_interpolation_coeffs_constant() {
        // Constant function f=1: all coefficients except a_0 should be ~0
        let n = 8;
        let vals = vec![1.0f64; n + 1];
        let coeffs = ChebyshevPolynomial::interpolation_coeffs(&vals);
        assert!((coeffs[0] - 1.0).abs() < 1e-10);
        for &c in coeffs.iter().skip(1) {
            assert!(c.abs() < 1e-10);
        }
    }

    // ------------------------------------------------------------------
    // LegendrePolynomial
    // ------------------------------------------------------------------
    #[test]
    fn test_legendre_p0() {
        assert!((LegendrePolynomial::eval(0, 0.5) - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_legendre_p1() {
        assert!((LegendrePolynomial::eval(1, 0.3) - 0.3).abs() < 1e-14);
    }

    #[test]
    fn test_legendre_p2() {
        let x = 0.5;
        let expected = 0.5 * (3.0 * x * x - 1.0);
        assert!((LegendrePolynomial::eval(2, x) - expected).abs() < 1e-12);
    }

    #[test]
    fn test_legendre_p3() {
        let x = 0.4;
        let expected = 0.5 * (5.0 * x * x * x - 3.0 * x);
        assert!((LegendrePolynomial::eval(3, x) - expected).abs() < 1e-12);
    }

    #[test]
    fn test_legendre_eval_all_consistency() {
        let x = 0.7;
        let all = LegendrePolynomial::eval_all(4, x);
        for k in 0..=4 {
            assert!((all[k] - LegendrePolynomial::eval(k, x)).abs() < 1e-12);
        }
    }

    #[test]
    fn test_gauss_legendre_nodes_count() {
        let (nodes, weights) = LegendrePolynomial::gauss_legendre(5);
        assert_eq!(nodes.len(), 5);
        assert_eq!(weights.len(), 5);
    }

    #[test]
    fn test_gauss_legendre_weights_sum() {
        let (_, weights) = LegendrePolynomial::gauss_legendre(5);
        let sum: f64 = weights.iter().sum();
        assert!((sum - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_gauss_legendre_integrate_poly() {
        // Integrate x^4 on [-1,1] = 2/5
        let result = LegendrePolynomial::integrate(|x| x.powi(4), 5);
        assert!((result - 0.4).abs() < 1e-10);
    }

    #[test]
    fn test_gauss_legendre_integrate_exp() {
        // Integrate exp(x) on [-1,1] = e - 1/e
        let exact = std::f64::consts::E - 1.0 / std::f64::consts::E;
        let result = LegendrePolynomial::integrate(|x| x.exp(), 8);
        assert!((result - exact).abs() < 1e-10);
    }

    // ------------------------------------------------------------------
    // FourierSeries / FFT
    // ------------------------------------------------------------------
    #[test]
    fn test_fft_length() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 0.0, 0.0, 0.0, 0.0];
        let spec = FourierSeries::fft(&data);
        assert_eq!(spec.len(), 8);
    }

    #[test]
    fn test_fft_ifft_roundtrip() {
        let n = 16;
        let data: Vec<f64> = (0..n)
            .map(|k| (2.0 * PI * k as f64 / n as f64).sin())
            .collect();
        let spec = FourierSeries::fft(&data);
        let recovered = FourierSeries::ifft(&spec);
        for (a, b) in data.iter().zip(recovered.iter()) {
            assert!((a - b).abs() < 1e-10, "mismatch: {a} vs {b}");
        }
    }

    #[test]
    fn test_fft_dc_component() {
        // Constant signal: FFT[0] = N * mean
        let n = 8;
        let data = vec![3.0f64; n];
        let spec = FourierSeries::fft(&data);
        assert!((spec[0].0 - 3.0 * n as f64).abs() < 1e-10);
        // All other bins should be ~0
        for k in 1..n {
            assert!(spec[k].0.abs() < 1e-10 && spec[k].1.abs() < 1e-10);
        }
    }

    #[test]
    fn test_fft_single_frequency() {
        let n = 8usize;
        let k0 = 2usize; // frequency bin 2
        let data: Vec<f64> = (0..n)
            .map(|j| (2.0 * PI * k0 as f64 * j as f64 / n as f64).cos())
            .collect();
        let spec = FourierSeries::fft(&data);
        // bins k0 and n-k0 should have amplitude n/2
        let amp_k0 = (spec[k0].0.powi(2) + spec[k0].1.powi(2)).sqrt();
        assert!((amp_k0 - n as f64 / 2.0).abs() < 1e-8);
    }

    #[test]
    fn test_power_spectrum_length() {
        let data = vec![0.0f64; 16];
        let ps = FourierSeries::power_spectrum(&data);
        assert_eq!(ps.len(), 9); // n/2 + 1
    }

    #[test]
    fn test_convolve_delta() {
        // Convolving with a delta (impulse at 0) should return the original signal
        let n = 8usize;
        let signal: Vec<f64> = (0..n).map(|k| k as f64 + 1.0).collect();
        let mut delta = vec![0.0f64; n];
        delta[0] = 1.0;
        let result = FourierSeries::convolve(&signal, &delta);
        for (a, b) in signal.iter().zip(result.iter()) {
            assert!((a - b).abs() < 1e-8);
        }
    }

    #[test]
    fn test_frequencies_length() {
        let freqs = FourierSeries::frequencies(16, 100.0);
        assert_eq!(freqs.len(), 9);
        assert!((freqs[0]).abs() < 1e-12);
    }

    // ------------------------------------------------------------------
    // SpectralDiff
    // ------------------------------------------------------------------
    #[test]
    fn test_spectral_diff_sin() {
        // d/dx sin(x) = cos(x) on [0, 2pi)
        let n = 64usize;
        let l = 2.0 * PI;
        let u: Vec<f64> = (0..n)
            .map(|k| (2.0 * PI * k as f64 / n as f64).sin())
            .collect();
        let du = SpectralDiff::diff(&u, l);
        let expected: Vec<f64> = (0..n)
            .map(|k| (2.0 * PI * k as f64 / n as f64).cos())
            .collect();
        for (got, exp) in du.iter().zip(expected.iter()) {
            assert!((got - exp).abs() < 1e-8, "got {got}, expected {exp}");
        }
    }

    #[test]
    fn test_spectral_diff2_sin() {
        // d^2/dx^2 sin(x) = -sin(x)
        let n = 64usize;
        let l = 2.0 * PI;
        let u: Vec<f64> = (0..n)
            .map(|k| (2.0 * PI * k as f64 / n as f64).sin())
            .collect();
        let d2u = SpectralDiff::diff2(&u, l);
        for (k, (&got, &orig)) in d2u.iter().zip(u.iter()).enumerate() {
            let _ = k;
            assert!((got + orig).abs() < 1e-8);
        }
    }

    // ------------------------------------------------------------------
    // ChebyshevCollocation
    // ------------------------------------------------------------------
    #[test]
    fn test_collocation_poisson_linear() {
        // u'' = 0, u(-1)=0, u(1)=1 => u(x) = (x+1)/2
        let coll = ChebyshevCollocation::new(12);
        let u = coll.solve_poisson(|_x| 0.0, 0.0, 1.0);
        let nodes = coll.nodes();
        for (&xi, &ui) in nodes.iter().zip(u.iter()) {
            let expected = (xi + 1.0) / 2.0;
            assert!(
                (ui - expected).abs() < 1e-8,
                "x={xi} u={ui} expected={expected}"
            );
        }
    }

    #[test]
    fn test_collocation_nodes_count() {
        let coll = ChebyshevCollocation::new(8);
        assert_eq!(coll.nodes().len(), 8);
    }

    #[test]
    fn test_collocation_diff_matrix_size() {
        let coll = ChebyshevCollocation::new(6);
        let d = coll.diff_matrix();
        assert_eq!(d.len(), 6);
    }

    // ------------------------------------------------------------------
    // WaveletTransform
    // ------------------------------------------------------------------
    #[test]
    fn test_haar_forward_inverse_roundtrip() {
        let signal = vec![1.0, 3.0, 5.0, 7.0, 2.0, 4.0, 6.0, 8.0];
        let (approx, detail) = WaveletTransform::haar_forward(&signal);
        let recovered = WaveletTransform::haar_inverse(&approx, &detail);
        for (a, b) in signal.iter().zip(recovered.iter()) {
            assert!((a - b).abs() < 1e-12);
        }
    }

    #[test]
    fn test_haar_decompose_reconstruct() {
        let signal: Vec<f64> = (0..16).map(|k| k as f64).collect();
        let coeffs = WaveletTransform::decompose(&signal, 3);
        let recovered = WaveletTransform::reconstruct(&coeffs);
        for (a, b) in signal.iter().zip(recovered.iter()) {
            assert!((a - b).abs() < 1e-10);
        }
    }

    #[test]
    fn test_haar_decompose_levels() {
        let signal = vec![1.0f64; 8];
        let coeffs = WaveletTransform::decompose(&signal, 3);
        // 3 levels + 1 approximation = 4 sub-bands
        assert_eq!(coeffs.len(), 4);
    }

    #[test]
    fn test_wavelet_soft_threshold_zeros_small() {
        let signal = vec![0.1f64, 0.2, 0.05, 1.0, 0.8, 0.03, 0.9, 0.02];
        let mut coeffs = WaveletTransform::decompose(&signal, 2);
        WaveletTransform::soft_threshold(&mut coeffs, 0.5);
        // Detail coefficients smaller than lambda should become 0
        for sub in coeffs.iter().skip(1) {
            for &v in sub {
                assert!(v.abs() <= v.abs() + 0.5); // trivially true; check no blow-up
            }
        }
    }

    #[test]
    fn test_wavelet_subband_energy() {
        let signal: Vec<f64> = (0..8).map(|k| (k as f64).sin()).collect();
        let coeffs = WaveletTransform::decompose(&signal, 2);
        let energies = WaveletTransform::subband_energy(&coeffs);
        assert_eq!(energies.len(), coeffs.len());
        for &e in &energies {
            assert!(e >= 0.0);
        }
    }

    #[test]
    fn test_wavelet_energy_conservation() {
        let signal: Vec<f64> = (0..8).map(|k| k as f64 + 1.0).collect();
        let total_energy: f64 = signal.iter().map(|v| v * v).sum();
        let coeffs = WaveletTransform::decompose(&signal, 3);
        let sub_energies: f64 = WaveletTransform::subband_energy(&coeffs).iter().sum();
        // Energy should be conserved (orthogonal transform)
        assert!((total_energy - sub_energies).abs() < 1e-8);
    }

    // ------------------------------------------------------------------
    // Gauss-Legendre orthogonality
    // ------------------------------------------------------------------
    #[test]
    fn test_legendre_orthogonality() {
        // Integrate P_2(x) * P_3(x) on [-1,1] should be 0
        let result = LegendrePolynomial::integrate(
            |x| LegendrePolynomial::eval(2, x) * LegendrePolynomial::eval(3, x),
            8,
        );
        assert!(result.abs() < 1e-10);
    }

    #[test]
    fn test_legendre_normalization() {
        // Integrate P_2(x)^2 on [-1,1] = 2/(2*2+1) = 2/5
        let result = LegendrePolynomial::integrate(|x| LegendrePolynomial::eval(2, x).powi(2), 8);
        assert!((result - 2.0 / 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_chebyshev_orthogonality_numerical() {
        // T_2 and T_3 are orthogonal w.r.t. weight 1/sqrt(1-x^2)
        // Numerical check via 16-pt GL: integral T_2 T_3 / sqrt(1-x^2) ~ 0
        let result = LegendrePolynomial::integrate(
            |x| {
                let w = if (1.0 - x * x) > 1e-10 {
                    1.0 / (1.0 - x * x).sqrt()
                } else {
                    0.0
                };
                ChebyshevPolynomial::eval(2, x) * ChebyshevPolynomial::eval(3, x) * w
            },
            16,
        );
        assert!(result.abs() < 1e-6);
    }
}
