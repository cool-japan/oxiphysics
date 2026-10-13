// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Fractional calculus module for the OxiPhysics engine.
//!
//! Implements key fractional calculus operators and models used in anomalous
//! diffusion, viscoelasticity, and non-local dynamics:
//!
//! - **Riemann-Liouville** fractional integral and derivative
//! - **Caputo** fractional derivative (initial-value-friendly form)
//! - **Grünwald-Letnikov** discrete approximation
//! - **Mittag-Leffler** function E_{α,β}(z) (one- and two-parameter)
//! - **Fractional differential equations** (FDE) via predictor-corrector
//! - **Anomalous diffusion**: subdiffusion (α<1) and superdiffusion (α>1)
//! - **Riesz fractional derivative** (symmetric, space-fractional)
//! - **Fractional Laplacian** approximation on uniform grids
//! - **Memory kernels**: power-law, exponential, Mittag-Leffler
//! - **Fractional oscillator** (Bagley-Torvik, fractional spring-dashpot)
//! - **Lévy stable distributions**: characteristic function and sampling

use crate::error::PhysicsError;
use std::f64::consts::PI;

// ─── Internal LCG RNG ────────────────────────────────────────────────────────

/// Minimal linear congruential generator for reproducible Monte Carlo.
struct Lcg {
    state: u64,
}

impl Lcg {
    fn new(seed: u64) -> Self {
        Self {
            state: seed.wrapping_add(1),
        }
    }
    fn next_u64(&mut self) -> u64 {
        self.state = self
            .state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.state
    }
    fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
    fn normal(&mut self) -> f64 {
        let u1 = self.uniform().max(1e-300);
        let u2 = self.uniform();
        (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos()
    }
}

// ============================================================================
// 1. Gamma Function
// ============================================================================

/// Lanczos approximation of the Gamma function Γ(x) for x > 0.
///
/// Accurate to ~15 significant digits for real positive arguments.
pub fn gamma(x: f64) -> f64 {
    if x <= 0.0 {
        // Reflection formula: Γ(x)Γ(1-x) = π/sin(πx)
        let sinpix = (PI * x).sin();
        if sinpix.abs() < 1e-300 {
            return f64::INFINITY;
        }
        return PI / (sinpix * gamma(1.0 - x));
    }
    if x < 0.5 {
        return PI / ((PI * x).sin() * gamma(1.0 - x));
    }
    // Lanczos coefficients g=7
    let g = 7.0_f64;
    let c = [
        0.999_999_999_999_809_9_f64,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    let z = x - 1.0;
    let mut sum = c[0];
    for (i, &ci) in c[1..].iter().enumerate() {
        sum += ci / (z + (i as f64) + 1.0);
    }
    let t = z + g + 0.5;
    (2.0 * PI).sqrt() * t.powf(z + 0.5) * (-t).exp() * sum
}

/// Log-Gamma: ln Γ(x) for x > 0.
///
/// Numerically stable version using the Lanczos result.
pub fn log_gamma(x: f64) -> f64 {
    if x <= 0.0 {
        return f64::INFINITY;
    }
    gamma(x).abs().ln()
}

// ============================================================================
// 2. Mittag-Leffler Function
// ============================================================================

/// One-parameter Mittag-Leffler function E_α(z) = Σ z^k / Γ(αk+1).
///
/// Converges for |z| small; use at most `max_terms` terms.
/// For z < 0 and 0 < α ≤ 1 this is the standard anomalous relaxation kernel.
pub fn mittag_leffler_1(alpha: f64, z: f64, max_terms: usize) -> f64 {
    assert!(alpha > 0.0, "alpha must be positive");
    let mut sum = 0.0_f64;
    let mut z_pow = 1.0_f64; // z^k
    for k in 0..max_terms {
        let g = gamma(alpha * k as f64 + 1.0);
        if g.is_infinite() || g.is_nan() {
            break;
        }
        let term = z_pow / g;
        sum += term;
        if term.abs() < 1e-15 * sum.abs().max(1e-300) {
            break;
        }
        z_pow *= z;
    }
    sum
}

/// Two-parameter Mittag-Leffler function E_{α,β}(z) = Σ z^k / Γ(αk+β).
///
/// Generalises the one-parameter version: E_{α,1}(z) = E_α(z).
/// Used in the solution of fractional differential equations.
pub fn mittag_leffler_2(alpha: f64, beta: f64, z: f64, max_terms: usize) -> f64 {
    assert!(alpha > 0.0, "alpha must be positive");
    assert!(beta > 0.0, "beta must be positive");
    let mut sum = 0.0_f64;
    let mut z_pow = 1.0_f64;
    for k in 0..max_terms {
        let g = gamma(alpha * k as f64 + beta);
        if g.is_infinite() || g.is_nan() {
            break;
        }
        let term = z_pow / g;
        sum += term;
        if term.abs() < 1e-15 * sum.abs().max(1e-300) {
            break;
        }
        z_pow *= z;
    }
    sum
}

/// Derivative of E_{α,1}(z) with respect to z, computed by differentiating
/// the series term by term: d/dz E_α(z) = E_{α,α}(z).
///
/// Useful for FDE Jacobians and sensitivity analysis.
pub fn mittag_leffler_derivative(alpha: f64, z: f64, max_terms: usize) -> f64 {
    mittag_leffler_2(alpha, alpha, z, max_terms)
}

// ============================================================================
// 3. Riemann-Liouville Fractional Integral
// ============================================================================

/// Riemann-Liouville fractional integral I^α\[f\] at point `t_n` using the
/// product-integration rule on a uniform grid with step `h`.
///
/// I^α\[f\](t_n) = (1/Γ(α)) Σ_{j=0}^{n} w_j * f(t_j)
///
/// where the weights are the standard RL quadrature coefficients.
///
/// # Arguments
/// * `f_values` – sampled function values f(t_0), …, f(t_n)
/// * `h`        – uniform step size
/// * `alpha`    – fractional order (0 < α < 1 for subdiffusion)
pub fn rl_integral(f_values: &[f64], h: f64, alpha: f64) -> f64 {
    assert!(alpha > 0.0, "alpha must be positive");
    let n = f_values.len();
    if n == 0 {
        return 0.0;
    }
    // Piecewise constant quadrature: I^α[f](t_n) ≈ (h^α/Γ(α+1)) Σ b_j f_j
    // where b_j = (n-j)^α - (n-j-1)^α  (rectangular rule, Γ(α+1) denominator)
    let g_alpha1 = gamma(alpha + 1.0);
    let mut sum = 0.0_f64;
    let n_idx = n - 1;
    for (j, &fv) in f_values.iter().enumerate() {
        let k = n_idx - j; // k = n-1-j goes from n-1 down to 0
        let b = (k as f64 + 1.0).powf(alpha) - (k as f64).powf(alpha);
        sum += b * fv;
    }
    h.powf(alpha) * sum / g_alpha1
}

/// Riemann-Liouville fractional integral I^α\[f\] evaluated at all grid points.
///
/// Returns a `Vec`f64` of length equal to `f_values`.
pub fn rl_integral_series(f_values: &[f64], h: f64, alpha: f64) -> Vec<f64> {
    (1..=f_values.len())
        .map(|n| rl_integral(&f_values[..n], h, alpha))
        .collect()
}

// ============================================================================
// 4. Riemann-Liouville Fractional Derivative
// ============================================================================

/// Grünwald-Letnikov coefficients g_j^(α) = (-1)^j C(α,j).
///
/// Used in both GL and RL discrete approximations.
fn gl_coeff(alpha: f64, j: usize) -> f64 {
    // Recurrence: g_0=1, g_j = g_{j-1} * (j - 1 - α) / j
    let mut g = 1.0_f64;
    for k in 1..=j {
        g *= (k as f64 - 1.0 - alpha) / k as f64;
    }
    g
}

/// Riemann-Liouville fractional derivative D^α\[f\] at t_n via the
/// Grünwald-Letnikov shifted approximation (first-order accurate).
///
/// D^α\[f\](t_n) ≈ h^{-α} Σ_{j=0}^{n} g_j * f(t_{n-j})
///
/// # Arguments
/// * `f_values` – f(t_0), …, f(t_n)
/// * `h`        – uniform time step
/// * `alpha`    – fractional order 0 < α < 2
pub fn rl_derivative(f_values: &[f64], h: f64, alpha: f64) -> f64 {
    assert!(alpha > 0.0 && alpha < 2.0, "alpha must be in (0,2)");
    let n = f_values.len();
    if n == 0 {
        return 0.0;
    }
    let mut sum = 0.0_f64;
    for j in 0..n {
        let g = gl_coeff(alpha, j);
        sum += g * f_values[n - 1 - j];
    }
    h.powf(-alpha) * sum
}

/// Riemann-Liouville fractional derivative at all grid points.
pub fn rl_derivative_series(f_values: &[f64], h: f64, alpha: f64) -> Vec<f64> {
    (1..=f_values.len())
        .map(|n| rl_derivative(&f_values[..n], h, alpha))
        .collect()
}

// ============================================================================
// 5. Caputo Fractional Derivative
// ============================================================================

/// Caputo fractional derivative ^C D^α\[f\](t_n) for 0 < α < 1.
///
/// The Caputo definition is preferable for initial-value problems because it
/// uses classical integer-order initial conditions.
///
/// ^C D^α\[f\](t_n) = (1/Γ(1-α)) Σ_{j=0}^{n-1} (f(t_{j+1})-f(t_j))/h * b_{n-1-j}
///
/// where b_k = ((k+1)^{1-α} - k^{1-α}) / (1-α) ... simplified via the L1 scheme.
///
/// # Arguments
/// * `f_values` – f(t_0), …, f(t_n)
/// * `h`        – uniform time step
/// * `alpha`    – fractional order 0 < α < 1
pub fn caputo_derivative(f_values: &[f64], h: f64, alpha: f64) -> f64 {
    assert!(alpha > 0.0 && alpha < 1.0, "alpha must be in (0,1)");
    let n = f_values.len();
    if n < 2 {
        return 0.0;
    }
    let mu = 1.0 - alpha;
    let g = gamma(mu + 1.0); // = Γ(2-α)
    let mut sum = 0.0_f64;
    let m = n - 1;
    for j in 0..m {
        let df = f_values[j + 1] - f_values[j];
        let b = ((m - j) as f64).powf(mu) - ((m - j - 1) as f64).powf(mu);
        sum += b * df;
    }
    sum / (h.powf(alpha) * g)
}

/// Caputo fractional derivative computed at every grid point.
pub fn caputo_derivative_series(f_values: &[f64], h: f64, alpha: f64) -> Vec<f64> {
    (2..=f_values.len())
        .map(|n| caputo_derivative(&f_values[..n], h, alpha))
        .collect()
}

// ============================================================================
// 6. Grünwald-Letnikov Definition
// ============================================================================

/// Grünwald-Letnikov fractional derivative of order α at t_n.
///
/// Direct definition: D^α\[f\](t_n) = lim_{h→0} h^{-α} Σ_{j=0}^{n} g_j f_{n-j}
///
/// Functionally identical to `rl_derivative` on a uniform grid; provided as a
/// named implementation reflecting the GL definition explicitly.
pub fn gl_derivative(f_values: &[f64], h: f64, alpha: f64) -> f64 {
    rl_derivative(f_values, h, alpha)
}

/// All GL coefficients g_0 … g_n for order α.
///
/// Efficient computation via the recurrence relation.
pub fn gl_coefficients(alpha: f64, n: usize) -> Vec<f64> {
    let mut coeffs = Vec::with_capacity(n + 1);
    let mut g = 1.0_f64;
    coeffs.push(g);
    for k in 1..=n {
        g *= (k as f64 - 1.0 - alpha) / k as f64;
        coeffs.push(g);
    }
    coeffs
}

// ============================================================================
// 7. Riesz Fractional Derivative
// ============================================================================

/// Riesz fractional derivative ∂^α/∂|x|^α on a uniform 1-D grid.
///
/// The Riesz derivative is symmetric and non-local.  It is the sum of the
/// left- and right-sided Grünwald-Letnikov derivatives,
///
/// (∂^α f / ∂|x|^α)_j ≈ c_α h^{-α} \[ Σ_{k=0}^{j} g_k f_{j-k} + Σ_{k=0}^{n-1-j} g_k f_{j+k} \]
///
/// where c_α = -1/(2 cos(πα/2)) and g_k = (-1)^k C(α, k) are the GL
/// coefficients.  The left-sided sum runs over every grid point to the left of
/// `x_j` (including itself) and the right-sided sum over every grid point to its
/// right, so symmetric data produce symmetric output.  The scheme is first-order
/// accurate in `h`.  Requires 0 < α < 2, α ≠ 1 (the factor c_α is singular at
/// α = 1; for α = 2s ≤ 2 including α = 1 use [`fractional_laplacian_1d`], which
/// is second-order accurate).
///
/// Boundary values outside the domain are treated as zero (zero-padding).
///
/// # Arguments
/// * `f`     – spatial function values on uniform grid
/// * `h`     – spatial step
/// * `alpha` – fractional order 0 < α < 2, α ≠ 1
pub fn riesz_derivative(f: &[f64], h: f64, alpha: f64) -> Vec<f64> {
    assert!(
        alpha > 0.0 && alpha < 2.0 && (alpha - 1.0).abs() > 1e-12,
        "alpha must be in (0,2) and not equal to 1"
    );
    let n = f.len();
    if n == 0 {
        return Vec::new();
    }
    let c = -1.0 / (2.0 * (PI * alpha / 2.0).cos());
    // g_0 … g_{n-1}: the longest sum (k = n-1) is the far end of the grid.
    let coeffs = gl_coefficients(alpha, n - 1);
    let scale = c / h.powf(alpha);
    (0..n)
        .map(|j| {
            // Left-sided: Σ_{k=0}^{j} g_k f_{j-k}.
            let left: f64 = coeffs[..=j]
                .iter()
                .zip(f[..=j].iter().rev())
                .map(|(g, fv)| g * fv)
                .sum();
            // Right-sided: Σ_{k=0}^{n-1-j} g_k f_{j+k}.
            let right: f64 = coeffs[..n - j]
                .iter()
                .zip(f[j..].iter())
                .map(|(g, fv)| g * fv)
                .sum();
            scale * (left + right)
        })
        .collect()
}

// ============================================================================
// 8. Fractional Laplacian
// ============================================================================

/// Coefficients `g_0 … g_{n-1}` of the centred (Ortigueira) fractional
/// difference of order `alpha = 2s`.
///
/// `g_k = (-1)^k Γ(α+1) / (Γ(α/2-k+1) Γ(α/2+k+1))`, computed through
/// `g_0 = Γ(α+1) / Γ(α/2+1)²` and `g_{k+1} = g_k (k - α/2) / (k + α/2 + 1)`.
/// The Fourier symbol of `Σ_k g_{|k|} e^{-ikθ}` is `(2 sin(θ/2))^α`, which is
/// `|θ|^α (1 + O(θ²))` and therefore a second-order accurate approximation of
/// `|ξ|^α`, the symbol of `(-Δ)^{α/2}`.
fn centred_difference_coefficients(alpha: f64, n: usize) -> Vec<f64> {
    let half = alpha / 2.0;
    let mut coeffs = Vec::with_capacity(n);
    if n == 0 {
        return coeffs;
    }
    let mut g = gamma(alpha + 1.0) / (gamma(half + 1.0) * gamma(half + 1.0));
    coeffs.push(g);
    for k in 1..n {
        g *= ((k - 1) as f64 - half) / (k as f64 + half);
        coeffs.push(g);
    }
    coeffs
}

/// Apply the symmetric Toeplitz operator `out_j = scale · Σ_i g_{|i-j|} line_i`
/// (values outside `line` are zero).  `g` must hold at least `line.len()`
/// coefficients.
fn apply_centred_difference(g: &[f64], line: &[f64], scale: f64) -> Vec<f64> {
    let n = line.len();
    (0..n)
        .map(|j| {
            // Σ_{k=1}^{j} g_k f_{j-k}
            let below: f64 = g[1..=j]
                .iter()
                .zip(line[..j].iter().rev())
                .map(|(gk, fv)| gk * fv)
                .sum();
            // Σ_{k=1}^{n-1-j} g_k f_{j+k}
            let above: f64 = g[1..n - j]
                .iter()
                .zip(line[j + 1..].iter())
                .map(|(gk, fv)| gk * fv)
                .sum();
            scale * (g[0] * line[j] + below + above)
        })
        .collect()
}

/// Validate the grid spacing `h` and fractional power `s ∈ (0, 1]`.
fn validate_laplacian_arguments(h: f64, s: f64) -> Result<(), PhysicsError> {
    if !(h.is_finite() && h > 0.0) {
        return Err(PhysicsError::InvalidInput {
            field: "h".to_string(),
            reason: format!("grid spacing must be positive and finite, got {h}"),
        });
    }
    if !(s > 0.0 && s <= 1.0) {
        return Err(PhysicsError::InvalidInput {
            field: "s".to_string(),
            reason: format!("fractional power must lie in (0, 1], got {s}"),
        });
    }
    Ok(())
}

/// Fractional Laplacian (-Δ)^{s} on a 1-D uniform grid, 0 < s ≤ 1.
///
/// Evaluated with the centred fractional difference of Ortigueira,
///
/// (-Δ)^s f_j ≈ h^{-2s} Σ_{k=-(n-1)}^{n-1} g_{|k|} f_{j-k},
///
/// with `g_0 = Γ(2s+1)/Γ(s+1)²`, `g_{k+1} = g_k (k - s)/(k + s + 1)` and `f`
/// taken as zero outside the grid.  The discrete symbol `(2 sin(θ/2))^{2s} h^{-2s}`
/// approximates `|ξ|^{2s}` to second order, the operator is symmetric positive
/// semi-definite, and it is valid on the whole range including `s = 1/2`
/// (where the Riesz Grünwald-Letnikov formula is singular) and `s = 1`, where it
/// reduces to the classical three-point stencil `(-f_{j-1} + 2 f_j - f_{j+1})/h²`.
///
/// The direct summation costs O(n²).
///
/// # Arguments
/// * `f`     – values on uniform grid of length n
/// * `h`     – grid spacing (positive, finite)
/// * `s`     – fractional power, 0 < s ≤ 1
///
/// # Errors
/// [`PhysicsError::InvalidInput`] if `h` is not positive and finite or `s` is
/// outside `(0, 1]`.
pub fn fractional_laplacian_1d(f: &[f64], h: f64, s: f64) -> Result<Vec<f64>, PhysicsError> {
    validate_laplacian_arguments(h, s)?;
    let alpha = 2.0 * s;
    let g = centred_difference_coefficients(alpha, f.len());
    Ok(apply_centred_difference(&g, f, h.powf(-alpha)))
}

/// Axis-wise fractional Laplacian on a 2-D uniform grid (n_x × n_y):
/// `(-∂ₓ²)^s + (-∂ᵧ²)^s`.
///
/// Row-major layout.  The 1-D operator of [`fractional_laplacian_1d`] is
/// applied along every grid row and every grid column and the results are
/// summed, with zero extension outside the grid.  For `s = 1` this is exactly the
/// five-point Laplacian `-Δ`.  For `s < 1` it is **not** the isotropic
/// fractional Laplacian `(-Δ)^s`: its Fourier symbol is `|k_x|^{2s} + |k_y|^{2s}`
/// rather than `(k_x² + k_y²)^s`, so it is rotation-dependent (an "anisotropic" /
/// directional fractional operator).
///
/// # Arguments
/// * `f`    – row-major 2-D data, length n_x * n_y
/// * `nx`   – number of x grid points
/// * `ny`   – number of y grid points
/// * `h`    – uniform grid spacing (same in x and y)
/// * `s`    – fractional power 0 < s ≤ 1
///
/// # Errors
/// [`PhysicsError::InvalidInput`] if `f.len() != nx * ny`, if `h` is not
/// positive and finite, or if `s` is outside `(0, 1]`.
pub fn fractional_laplacian_2d(
    f: &[f64],
    nx: usize,
    ny: usize,
    h: f64,
    s: f64,
) -> Result<Vec<f64>, PhysicsError> {
    validate_laplacian_arguments(h, s)?;
    let expected = nx
        .checked_mul(ny)
        .ok_or_else(|| PhysicsError::InvalidInput {
            field: "nx*ny".to_string(),
            reason: format!("grid size {nx} x {ny} overflows usize"),
        })?;
    if f.len() != expected {
        return Err(PhysicsError::InvalidInput {
            field: "f".to_string(),
            reason: format!(
                "expected {expected} values for a {nx} x {ny} grid, got {}",
                f.len()
            ),
        });
    }
    let alpha = 2.0 * s;
    let scale = h.powf(-alpha);
    let gx = centred_difference_coefficients(alpha, nx);
    let gy = centred_difference_coefficients(alpha, ny);
    let mut result = vec![0.0_f64; expected];
    // x-direction: contiguous rows.
    if nx > 0 {
        for (row_in, row_out) in f.chunks_exact(nx).zip(result.chunks_exact_mut(nx)) {
            let d = apply_centred_difference(&gx, row_in, scale);
            for (out, dv) in row_out.iter_mut().zip(d.iter()) {
                *out += dv;
            }
        }
    }
    // y-direction: gather each column.
    for ix in 0..nx {
        let col: Vec<f64> = f.iter().skip(ix).step_by(nx).copied().collect();
        let d = apply_centred_difference(&gy, &col, scale);
        for (iy, dv) in d.iter().enumerate() {
            result[iy * nx + ix] += dv;
        }
    }
    Ok(result)
}

// ============================================================================
// 9. Memory Kernels
// ============================================================================

/// Power-law memory kernel K(t) = t^{α-1} / Γ(α).
///
/// This is the kernel of the Riemann-Liouville fractional integral.
/// For 0 < α < 1 it is weakly singular at t = 0.
pub fn memory_kernel_power_law(t: f64, alpha: f64) -> f64 {
    assert!(alpha > 0.0);
    if t <= 0.0 {
        return 0.0;
    }
    t.powf(alpha - 1.0) / gamma(alpha)
}

/// Exponential memory kernel K(t) = λ e^{-λt}.
///
/// Corresponds to standard Markovian (α → 1) limit of fractional dynamics.
pub fn memory_kernel_exponential(t: f64, lambda: f64) -> f64 {
    if t < 0.0 {
        return 0.0;
    }
    lambda * (-lambda * t).exp()
}

/// Mittag-Leffler memory kernel K(t) = t^{α-1} E_{α,α}(-λ t^α).
///
/// Appears in the solution of fractional relaxation equations.
/// Interpolates between exponential (α=1) and power-law decay.
pub fn memory_kernel_mittag_leffler(t: f64, alpha: f64, lambda: f64) -> f64 {
    assert!(alpha > 0.0 && alpha <= 1.0);
    if t <= 0.0 {
        return 0.0;
    }
    t.powf(alpha - 1.0) * mittag_leffler_2(alpha, alpha, -lambda * t.powf(alpha), 100)
}

/// Compute the convolution integral of a memory kernel K with function f on
/// a uniform grid [0, (n-1)*h] using the trapezoidal rule.
///
/// Returns the convolution value at the last grid point t_n.
///
/// (K * f)(t_n) = ∫_0^{t_n} K(t_n - s) f(s) ds
pub fn memory_convolution(kernel_vals: &[f64], f_vals: &[f64], h: f64) -> f64 {
    assert_eq!(kernel_vals.len(), f_vals.len());
    let n = kernel_vals.len();
    if n == 0 {
        return 0.0;
    }
    // Trapezoidal rule
    let mut sum = 0.5 * (kernel_vals[0] * f_vals[n - 1] + kernel_vals[n - 1] * f_vals[0]);
    for i in 1..n - 1 {
        sum += kernel_vals[i] * f_vals[n - 1 - i];
    }
    h * sum
}

// ============================================================================
// 10. Anomalous Diffusion
// ============================================================================

/// Parameters for an anomalous diffusion model.
///
/// Subdiffusion: 0 < α < 1, mean-squared displacement ⟨x²⟩ ~ t^α.
/// Normal diffusion: α = 1, ⟨x²⟩ ~ t.
/// Superdiffusion: 1 < α < 2, ⟨x²⟩ ~ t^α.
#[derive(Debug, Clone, Copy)]
pub struct AnomalousDiffusionParams {
    /// Fractional order α ∈ (0, 2).
    pub alpha: f64,
    /// Generalized diffusion coefficient K_α [m² s^{-α}].
    pub diffusivity: f64,
}

impl AnomalousDiffusionParams {
    /// Create new anomalous diffusion parameters.
    pub fn new(alpha: f64, diffusivity: f64) -> Self {
        assert!(alpha > 0.0 && alpha < 2.0, "alpha must be in (0,2)");
        assert!(diffusivity > 0.0, "diffusivity must be positive");
        Self { alpha, diffusivity }
    }

    /// Mean squared displacement at time t: ⟨x²⟩(t) = 2 K_α t^α / Γ(1+α).
    pub fn msd(&self, t: f64) -> f64 {
        if t <= 0.0 {
            return 0.0;
        }
        2.0 * self.diffusivity * t.powf(self.alpha) / gamma(1.0 + self.alpha)
    }

    /// Diffusion exponent (anomalous scaling exponent).
    pub fn anomalous_exponent(&self) -> f64 {
        self.alpha
    }

    /// Check if the process is subdiffusive.
    pub fn is_subdiffusion(&self) -> bool {
        self.alpha < 1.0
    }

    /// Check if the process is superdiffusive.
    pub fn is_superdiffusion(&self) -> bool {
        self.alpha > 1.0
    }
}

/// Subdiffusion probability density (Green's function) via Mittag-Leffler.
///
/// For the fractional diffusion equation with Caputo time derivative:
/// ∂^α C / ∂t^α = K_α ∂² C / ∂x²
///
/// The fundamental solution in an infinite domain is expressed through
/// the Fox H-function; this function returns the Gaussian approximation
/// valid for large |x| or small α deviations from 1.
///
/// # Arguments
/// * `x`      – spatial coordinate
/// * `t`      – time
/// * `params` – anomalous diffusion parameters
pub fn subdiffusion_pdf(x: f64, t: f64, params: &AnomalousDiffusionParams) -> f64 {
    if t <= 0.0 {
        return 0.0;
    }
    let sigma2 = params.msd(t);
    let _sigma = sigma2.sqrt().max(1e-300);
    // Gaussian approximation (exact for α=1, approximate otherwise)
    (-x * x / (2.0 * sigma2)).exp() / ((2.0 * PI * sigma2).sqrt())
}

// ============================================================================
// 11. Fractional Differential Equations (Predictor-Corrector)
// ============================================================================

/// Solve a scalar fractional initial-value problem using the Adams-Bashforth-
/// Moulton predictor-corrector method of Diethelm, Ford and Freed (2002).
///
/// Equation: ^C D^α y(t) = f(t, y(t)),  y(0) = y0,  0 < α ≤ 1.
///
/// With `h = t_end / n_steps`, `t_j = j h` and `f_j = f(t_j, y_j)`, the
/// equivalent Volterra equation `y(t) = y0 + (1/Γ(α)) ∫ (t-s)^{α-1} f(s, y(s)) ds`
/// is discretised by product integration.  The predictor is the fractional
/// Adams-Bashforth step
///
/// y^P_{n+1} = y0 + h^α/Γ(α+1) · Σ_{j=0}^{n} b_{n-j} f_j,
/// b_m = (m+1)^α - m^α,
///
/// and the corrector is the fractional Adams-Moulton step
///
/// y_{n+1} = y0 + h^α/Γ(α+2) · \[ f(t_{n+1}, y^P_{n+1}) + a_0(n) f_0 + Σ_{j=1}^{n} c_{n-j} f_j \],
/// a_0(n) = n^{α+1} - (n-α)(n+1)^α,
/// c_m = (m+2)^{α+1} + m^{α+1} - 2(m+1)^{α+1}.
///
/// For α = 1 this is the classical Euler-predictor / trapezoidal-corrector
/// (Heun) method.  The full history is used, so the cost is O(n_steps²) with two
/// evaluations of `f` per step (the value `f_j` is cached once).  The method
/// converges with order `min(2, 1 + α)` for smooth right-hand sides.
///
/// Returns `(t_values, y_values)` with `n_steps + 1` points.
///
/// # Arguments
/// * `f`       – right-hand side f(t, y)
/// * `y0`      – initial condition
/// * `t_end`   – final time
/// * `n_steps` – number of time steps
/// * `alpha`   – fractional order 0 < α ≤ 1
pub fn fde_predictor_corrector<F>(
    f: F,
    y0: f64,
    t_end: f64,
    n_steps: usize,
    alpha: f64,
) -> (Vec<f64>, Vec<f64>)
where
    F: Fn(f64, f64) -> f64,
{
    assert!(alpha > 0.0 && alpha <= 1.0);
    assert!(n_steps >= 1);
    let h = t_end / n_steps as f64;
    let t: Vec<f64> = (0..=n_steps).map(|k| k as f64 * h).collect();
    let mut y = vec![0.0_f64; n_steps + 1];
    y[0] = y0;
    // Cached right-hand-side values f_j = f(t_j, y_j).
    let mut fvals = vec![0.0_f64; n_steps + 1];
    fvals[0] = f(t[0], y0);

    // Weights that depend only on the lag m = n - j.
    let predictor_weights: Vec<f64> = (0..n_steps)
        .map(|m| (m as f64 + 1.0).powf(alpha) - (m as f64).powf(alpha))
        .collect();
    let corrector_weights: Vec<f64> = (0..n_steps)
        .map(|m| {
            let m = m as f64;
            (m + 2.0).powf(alpha + 1.0) + m.powf(alpha + 1.0) - 2.0 * (m + 1.0).powf(alpha + 1.0)
        })
        .collect();
    let h_alpha = h.powf(alpha);
    let predictor_scale = h_alpha / gamma(alpha + 1.0);
    let corrector_scale = h_alpha / gamma(alpha + 2.0);

    for n in 0..n_steps {
        // --- Predictor (fractional Adams-Bashforth) ---
        // Σ_{j=0}^{n} b_{n-j} f_j: pair f_0 … f_n with b_n … b_0.
        let predictor_sum: f64 = fvals[..=n]
            .iter()
            .zip(predictor_weights[..=n].iter().rev())
            .map(|(fj, b)| fj * b)
            .sum();
        let y_predicted = y0 + predictor_scale * predictor_sum;

        // --- Corrector (fractional Adams-Moulton) ---
        let nf = n as f64;
        let first_weight = nf.powf(alpha + 1.0) - (nf - alpha) * (nf + 1.0).powf(alpha);
        // Σ_{j=1}^{n} c_{n-j} f_j: pair f_1 … f_n with c_{n-1} … c_0.
        let history_sum: f64 = fvals[1..=n]
            .iter()
            .zip(corrector_weights[..n].iter().rev())
            .map(|(fj, c)| fj * c)
            .sum();
        let corrector_sum = f(t[n + 1], y_predicted) + first_weight * fvals[0] + history_sum;
        y[n + 1] = y0 + corrector_scale * corrector_sum;
        fvals[n + 1] = f(t[n + 1], y[n + 1]);
    }
    (t, y)
}

// ============================================================================
// 12. Fractional Oscillator
// ============================================================================

/// Fractional oscillator model: m D^{2α} x + b D^α x + k x = F(t).
///
/// For α = 1 this reduces to the classical damped harmonic oscillator.
/// For 0 < α < 1 the oscillator exhibits anomalous (sub-diffusive) behavior.
#[derive(Debug, Clone, Copy)]
pub struct FractionalOscillator {
    /// Mass m.
    pub mass: f64,
    /// Fractional damping coefficient b.
    pub damping: f64,
    /// Spring stiffness k.
    pub stiffness: f64,
    /// Fractional order α ∈ (0, 1].
    pub alpha: f64,
}

impl FractionalOscillator {
    /// Create a new fractional oscillator.
    pub fn new(mass: f64, damping: f64, stiffness: f64, alpha: f64) -> Self {
        assert!(alpha > 0.0 && alpha <= 1.0);
        Self {
            mass,
            damping,
            stiffness,
            alpha,
        }
    }

    /// Natural frequency ω_n for the integer-order limit (α=1).
    pub fn natural_frequency(&self) -> f64 {
        (self.stiffness / self.mass).sqrt()
    }

    /// Damping ratio ζ for the integer-order limit.
    pub fn damping_ratio(&self) -> f64 {
        self.damping / (2.0 * (self.mass * self.stiffness).sqrt())
    }

    /// Relaxation time τ = (m/k)^{1/(2α)}.
    pub fn relaxation_time(&self) -> f64 {
        (self.mass / self.stiffness).powf(1.0 / (2.0 * self.alpha))
    }

    /// Free response via Mittag-Leffler: x(t) = x0 E_{2α,1}(-ω² t^{2α}).
    ///
    /// Valid for underdamped (ζ < 1) or undamped (b=0) oscillators.
    pub fn free_response(&self, x0: f64, t: f64) -> f64 {
        if t <= 0.0 {
            return x0;
        }
        let omega2 = self.stiffness / self.mass;
        let two_alpha = 2.0 * self.alpha;
        x0 * mittag_leffler_1(two_alpha, -omega2 * t.powf(two_alpha), 80)
    }
}

/// Bagley-Torvik equation solver: m x'' + A D^{3/2} x + k x = F(t).
///
/// Returns time and displacement arrays for the Bagley-Torvik fractional model.
/// Uses a simplified Euler-fractional step as a first-order approximation.
///
/// # Arguments
/// * `mass`      – mass m
/// * `a_coeff`   – fractional damping coefficient A
/// * `stiffness` – spring stiffness k
/// * `force`     – external force function F(t)
/// * `x0`        – initial displacement
/// * `v0`        – initial velocity
/// * `t_end`     – final time
/// * `n_steps`   – number of steps
pub fn bagley_torvik<F>(
    mass: f64,
    a_coeff: f64,
    stiffness: f64,
    force: F,
    x0: f64,
    v0: f64,
    t_end: f64,
    n_steps: usize,
) -> (Vec<f64>, Vec<f64>)
where
    F: Fn(f64) -> f64,
{
    let h = t_end / n_steps as f64;
    let mut t = vec![0.0; n_steps + 1];
    let mut x = vec![0.0; n_steps + 1];
    let mut v = vec![0.0; n_steps + 1];
    x[0] = x0;
    v[0] = v0;
    for k in 0..n_steps {
        t[k + 1] = (k + 1) as f64 * h;
    }
    // Fractional GL coefficients for order 3/2
    let alpha_bt = 1.5_f64;
    let gl = gl_coefficients(alpha_bt, n_steps + 1);
    for n_idx in 0..n_steps {
        // Fractional term: D^{3/2} x ≈ h^{-3/2} Σ g_j x_{n-j}
        let mut frac_sum = 0.0_f64;
        for j in 0..=n_idx {
            frac_sum += gl[j] * x[n_idx - j];
        }
        let frac_deriv = h.powf(-alpha_bt) * frac_sum;
        // Acceleration from equation of motion
        let acc = (force(t[n_idx]) - a_coeff * frac_deriv - stiffness * x[n_idx]) / mass;
        v[n_idx + 1] = v[n_idx] + h * acc;
        x[n_idx + 1] = x[n_idx] + h * v[n_idx];
    }
    (t, x)
}

// ============================================================================
// 13. Lévy Stable Distributions
// ============================================================================

/// Parameters of a Lévy stable distribution S(α, β, c, μ).
///
/// The characteristic function is:
/// φ(θ) = exp(i μ θ - |cθ|^α (1 - iβ sign(θ) tan(πα/2))) for α ≠ 1
#[derive(Debug, Clone, Copy)]
pub struct LevyStableParams {
    /// Stability index α ∈ (0, 2].
    pub alpha: f64,
    /// Skewness β ∈ [-1, 1].
    pub beta: f64,
    /// Scale c > 0.
    pub scale: f64,
    /// Location μ.
    pub location: f64,
}

impl LevyStableParams {
    /// Create Lévy stable parameters.
    pub fn new(alpha: f64, beta: f64, scale: f64, location: f64) -> Self {
        assert!(alpha > 0.0 && alpha <= 2.0);
        assert!((-1.0..=1.0).contains(&beta));
        assert!(scale > 0.0);
        Self {
            alpha,
            beta,
            scale,
            location,
        }
    }

    /// Symmetric Lévy stable (β=0, μ=0).
    pub fn symmetric(alpha: f64, scale: f64) -> Self {
        Self::new(alpha, 0.0, scale, 0.0)
    }

    /// Standard Cauchy distribution (α=1, β=0).
    pub fn cauchy(scale: f64, location: f64) -> Self {
        Self::new(1.0, 0.0, scale, location)
    }

    /// Standard Gaussian distribution (α=2, β=0).
    pub fn gaussian(sigma: f64) -> Self {
        // For S(2,0,c,0): variance = 2c², so c = sigma/sqrt(2)
        Self::new(2.0, 0.0, sigma / 2.0_f64.sqrt(), 0.0)
    }
}

/// Sample from a Lévy stable distribution using the Chambers-Mallows-Stuck algorithm.
///
/// Returns a single sample from S(α, β, c, μ).
///
/// # Arguments
/// * `params` – Lévy stable parameters
/// * `seed`   – random seed for reproducibility
pub fn levy_stable_sample(params: &LevyStableParams, seed: u64) -> f64 {
    let mut rng = Lcg::new(seed);
    let alpha = params.alpha;
    let beta = params.beta;

    // Special case: Gaussian α=2
    if (alpha - 2.0).abs() < 1e-12 {
        return params.location + params.scale * 2.0_f64.sqrt() * rng.normal();
    }
    // Special case: Cauchy α=1, β=0
    if (alpha - 1.0).abs() < 1e-12 && beta.abs() < 1e-12 {
        let u = (PI * (rng.uniform() - 0.5)).tan();
        return params.location + params.scale * u;
    }

    // Chambers-Mallows-Stuck (1976) algorithm
    let u = PI * (rng.uniform() - 0.5); // uniform on (-π/2, π/2)
    let w = -rng.uniform().max(1e-300).ln(); // exponential(1)

    let sample = if (alpha - 1.0).abs() > 1e-12 {
        let xi = (PI / 2.0) * beta * (PI * alpha / 2.0).tan().atan() / (PI / 2.0);
        // Simplified CMS without the exact ξ shift for non-unit α
        let b = (PI / 2.0 * alpha).tan();
        let s = (1.0 + beta * beta * b * b).powf(1.0 / (2.0 * alpha));
        let phi0 = (beta * b).atan() / alpha;
        s * ((alpha * (u + phi0)).sin() / u.cos().powf(1.0 / alpha))
            * ((u - alpha * (u + phi0)).cos() / w).powf((1.0 - alpha) / alpha)
            - xi
    } else {
        // α = 1
        (PI / 2.0 + beta * u).tan() - beta * ((PI / 2.0 + beta * u) * w / (PI / 2.0 * u.cos())).ln()
    };

    params.location + params.scale * sample
}

/// Log characteristic function of a symmetric Lévy stable distribution.
///
/// ln φ(θ) = -|cθ|^α for the symmetric case (β=0, μ=0).
pub fn levy_log_char_fn_symmetric(theta: f64, alpha: f64, scale: f64) -> f64 {
    -(scale * theta.abs()).powf(alpha)
}

/// Compute the PDF of a symmetric Lévy stable distribution by numerical Fourier inversion.
///
/// Uses trapezoidal integration over frequency domain [-θ_max, θ_max].
///
/// # Arguments
/// * `x`       – evaluation point
/// * `alpha`   – stability index
/// * `scale`   – scale parameter
/// * `n_quad`  – number of quadrature points (must be even)
/// * `theta_max` – truncation of frequency domain
pub fn levy_stable_pdf(x: f64, alpha: f64, scale: f64, n_quad: usize, theta_max: f64) -> f64 {
    // Inverse Fourier transform: p(x) = (1/2π) ∫ φ(θ) e^{-iθx} dθ
    // For symmetric φ: integrate over [-θ_max, θ_max] and divide by 2π.
    let dtheta = 2.0 * theta_max / n_quad as f64;
    let mut sum = 0.0_f64;
    for k in 0..=n_quad {
        let theta = -theta_max + k as f64 * dtheta;
        let log_cf = levy_log_char_fn_symmetric(theta, alpha, scale);
        let weight = if k == 0 || k == n_quad { 0.5 } else { 1.0 };
        sum += weight * (log_cf.exp()) * (theta * x).cos();
    }
    sum * dtheta / (2.0 * PI)
}

// ============================================================================
// 14. Fractional Relaxation
// ============================================================================

/// Fractional relaxation equation: D^α y + y/τ = 0, y(0) = y0.
///
/// Exact solution: y(t) = y0 E_α(-(t/τ)^α).
///
/// Models anomalous (non-exponential) relaxation in complex systems.
///
/// # Arguments
/// * `y0`    – initial value
/// * `t`     – time
/// * `alpha` – fractional order 0 < α ≤ 1
/// * `tau`   – relaxation time constant
pub fn fractional_relaxation(y0: f64, t: f64, alpha: f64, tau: f64) -> f64 {
    if t < 0.0 {
        return y0;
    }
    y0 * mittag_leffler_1(alpha, -(t / tau).powf(alpha), 150)
}

/// Time array for fractional relaxation to fall below threshold.
///
/// Scans a time grid [0, t_max] and returns the first time where
/// |y(t)/y0| < threshold.  Returns None if threshold is never reached.
pub fn relaxation_time_threshold(
    alpha: f64,
    tau: f64,
    threshold: f64,
    t_max: f64,
    n_pts: usize,
) -> Option<f64> {
    for k in 0..=n_pts {
        let t = k as f64 * t_max / n_pts as f64;
        let y = fractional_relaxation(1.0, t, alpha, tau);
        if y.abs() < threshold {
            return Some(t);
        }
    }
    None
}

// ============================================================================
// 15. Subdiffusion / Superdiffusion Simulation
// ============================================================================

/// Single-particle subdiffusion trajectory via subordination.
///
/// Generates a discrete-time trajectory of a subdiffusing particle
/// using the random-walk subordination method: operational time is drawn
/// from a stable subordinator with index α.
///
/// # Arguments
/// * `n_steps`   – number of steps
/// * `alpha`     – subdiffusion exponent 0 < α < 1
/// * `diffusivity` – generalized diffusion coefficient
/// * `seed`      – random seed
///
/// Returns `(x, y)` position arrays.
pub fn subdiffusion_trajectory(
    n_steps: usize,
    alpha: f64,
    diffusivity: f64,
    seed: u64,
) -> (Vec<f64>, Vec<f64>) {
    assert!(alpha > 0.0 && alpha < 1.0);
    let mut rng = Lcg::new(seed);
    let mut x = vec![0.0_f64; n_steps + 1];
    let mut y = vec![0.0_f64; n_steps + 1];
    for k in 0..n_steps {
        // Operational time increment from a one-sided stable distribution
        // Approximated by a Weibull-like draw: dt_op ~ Gamma(alpha+1)^{1/alpha}
        let u = rng.uniform().max(1e-300);
        let dt_op = (-u.ln()).powf(1.0 / alpha);
        let sigma = (2.0 * diffusivity * dt_op).sqrt();
        x[k + 1] = x[k] + sigma * rng.normal();
        y[k + 1] = y[k] + sigma * rng.normal();
    }
    (x, y)
}

/// Superdiffusion trajectory via Lévy flights.
///
/// Steps are drawn from a symmetric Lévy stable distribution with index α ∈ (1, 2).
/// The long-tailed step lengths produce superdiffusive dynamics.
///
/// # Arguments
/// * `n_steps` – number of Lévy flight steps
/// * `alpha`   – stability index 1 < α < 2
/// * `scale`   – scale parameter c
/// * `seed`    – random seed
pub fn superdiffusion_trajectory(
    n_steps: usize,
    alpha: f64,
    scale: f64,
    seed: u64,
) -> (Vec<f64>, Vec<f64>) {
    assert!(alpha > 1.0 && alpha < 2.0);
    let params = LevyStableParams::symmetric(alpha, scale);
    let mut x = vec![0.0_f64; n_steps + 1];
    let mut y = vec![0.0_f64; n_steps + 1];
    for k in 0..n_steps {
        let dx = levy_stable_sample(&params, seed.wrapping_add(k as u64 * 2));
        let dy = levy_stable_sample(&params, seed.wrapping_add(k as u64 * 2 + 1));
        x[k + 1] = x[k] + dx;
        y[k + 1] = y[k] + dy;
    }
    (x, y)
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gamma_integer_values() {
        // Γ(n) = (n-1)!
        assert!((gamma(1.0) - 1.0).abs() < 1e-10);
        assert!((gamma(2.0) - 1.0).abs() < 1e-10);
        assert!((gamma(3.0) - 2.0).abs() < 1e-10);
        assert!((gamma(4.0) - 6.0).abs() < 1e-10);
        assert!((gamma(5.0) - 24.0).abs() < 1e-8);
    }

    #[test]
    fn test_gamma_half() {
        // Γ(1/2) = sqrt(π)
        let g = gamma(0.5);
        assert!((g - PI.sqrt()).abs() < 1e-10, "Γ(1/2) = {:.6}", g);
    }

    #[test]
    fn test_mittag_leffler_1_at_zero() {
        // E_α(0) = 1 for any α
        for alpha in [0.3, 0.5, 0.7, 1.0, 1.5] {
            let ml = mittag_leffler_1(alpha, 0.0, 50);
            assert!(
                (ml - 1.0).abs() < 1e-12,
                "E_{alpha}(0) = {:.6} for alpha={}",
                ml,
                alpha
            );
        }
    }

    #[test]
    fn test_mittag_leffler_1_reduces_to_exp() {
        // E_1(z) = e^z
        for z in [-1.0, 0.0, 0.5, 1.0] {
            let ml = mittag_leffler_1(1.0, z, 100);
            let exp = z.exp();
            assert!(
                (ml - exp).abs() < 1e-8,
                "E_1({}) = {:.6}, exp({}) = {:.6}",
                z,
                ml,
                z,
                exp
            );
        }
    }

    #[test]
    fn test_mittag_leffler_2_reduces_to_exp() {
        // E_{1,1}(z) = e^z
        for z in [-0.5, 0.0, 0.5] {
            let ml = mittag_leffler_2(1.0, 1.0, z, 100);
            let exp = z.exp();
            assert!(
                (ml - exp).abs() < 1e-8,
                "E_{{1,1}}({}) = {:.6}, e^z={:.6}",
                z,
                ml,
                exp
            );
        }
    }

    #[test]
    fn test_gl_coefficients_first_few() {
        // For α=1: g_0=1, g_1=-1, g_2=0, ...
        let coeffs = gl_coefficients(1.0, 3);
        assert!((coeffs[0] - 1.0).abs() < 1e-12);
        assert!((coeffs[1] - (-1.0)).abs() < 1e-12);
        assert!(coeffs[2].abs() < 1e-12, "g_2 for alpha=1: {:.6}", coeffs[2]);
    }

    #[test]
    fn test_rl_integral_constant_function() {
        // I^α[1](t) = t^α / Γ(α+1)
        let n = 100;
        let t_end = 1.0;
        let h = t_end / n as f64;
        let f_vals: Vec<f64> = vec![1.0; n + 1];
        let alpha = 0.5;
        let result = rl_integral(&f_vals, h, alpha);
        let exact = t_end.powf(alpha) / gamma(alpha + 1.0);
        assert!(
            (result - exact).abs() < 0.05,
            "RL integral: {:.6} vs {:.6}",
            result,
            exact
        );
    }

    #[test]
    fn test_caputo_derivative_power_function() {
        // ^C D^α [t^β] = Γ(β+1)/Γ(β-α+1) * t^{β-α}
        let alpha = 0.5;
        let beta = 1.5;
        let n = 200;
        let t_end = 2.0;
        let h = t_end / n as f64;
        let f_vals: Vec<f64> = (0..=n).map(|k| (k as f64 * h).powf(beta)).collect();
        let result = caputo_derivative(&f_vals, h, alpha);
        let t_n = t_end;
        let exact = gamma(beta + 1.0) / gamma(beta - alpha + 1.0) * t_n.powf(beta - alpha);
        let rel_err = (result - exact).abs() / exact.abs();
        assert!(
            rel_err < 0.05,
            "Caputo D^0.5[t^1.5]: {:.6} vs {:.6}, err={:.4}",
            result,
            exact,
            rel_err
        );
    }

    #[test]
    fn test_riesz_derivative_shape() {
        // Riesz derivative of a Gaussian bump should be smooth with zero mean
        let n = 64;
        let h = 0.1;
        let center = n as f64 * h / 2.0;
        let f: Vec<f64> = (0..n)
            .map(|i| {
                let x = i as f64 * h - center;
                (-x * x / 0.5).exp()
            })
            .collect();
        let d = riesz_derivative(&f, h, 0.5);
        assert_eq!(d.len(), n);
        // Check the result is finite
        for &v in &d {
            assert!(v.is_finite(), "Riesz derivative contains non-finite value");
        }
    }

    #[test]
    fn test_fractional_laplacian_1d() {
        let n = 32;
        let h = 0.2;
        let f: Vec<f64> = (0..n)
            .map(|i| {
                let x = (i as f64 - n as f64 / 2.0) * h;
                (-x * x).exp()
            })
            .collect();
        let lap = fractional_laplacian_1d(&f, h, 0.5).expect("valid arguments");
        assert_eq!(lap.len(), n);
        for &v in &lap {
            assert!(v.is_finite());
        }
    }

    #[test]
    fn test_memory_kernel_power_law() {
        let alpha = 0.5;
        // K(1.0) = 1 / Γ(0.5) = 1 / sqrt(π)
        let k = memory_kernel_power_law(1.0, alpha);
        let expected = 1.0 / gamma(alpha);
        assert!((k - expected).abs() < 1e-10);
    }

    #[test]
    fn test_memory_kernel_exponential() {
        let lambda = 2.0;
        let t = 1.0;
        let k = memory_kernel_exponential(t, lambda);
        let expected = lambda * (-lambda * t).exp();
        assert!((k - expected).abs() < 1e-12);
        // At t<0 should be zero
        assert_eq!(memory_kernel_exponential(-1.0, lambda), 0.0);
    }

    #[test]
    fn test_anomalous_diffusion_msd() {
        let params = AnomalousDiffusionParams::new(0.5, 1.0);
        // MSD at t=0 should be 0
        assert_eq!(params.msd(0.0), 0.0);
        // MSD should grow with time
        assert!(params.msd(1.0) > params.msd(0.5));
        // Check formula: 2 K t^α / Γ(1+α)
        let t = 2.0_f64;
        let expected = 2.0 * 1.0 * t.powf(0.5) / gamma(1.5);
        assert!((params.msd(t) - expected).abs() < 1e-10);
    }

    #[test]
    fn test_anomalous_diffusion_classification() {
        let sub = AnomalousDiffusionParams::new(0.5, 1.0);
        let normal = AnomalousDiffusionParams::new(1.0, 1.0);
        let sup = AnomalousDiffusionParams::new(1.5, 1.0);
        assert!(sub.is_subdiffusion());
        assert!(!sub.is_superdiffusion());
        assert!(!normal.is_subdiffusion());
        assert!(!normal.is_superdiffusion());
        assert!(sup.is_superdiffusion());
    }

    #[test]
    fn test_fractional_relaxation_initial_value() {
        // At t=0, y = y0
        let y = fractional_relaxation(3.0, 0.0, 0.5, 1.0);
        assert!((y - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_fractional_relaxation_decay() {
        // y(t) should decay towards 0 for large t
        let y0 = 1.0;
        let alpha = 0.7;
        let tau = 1.0;
        let y_large = fractional_relaxation(y0, 10.0, alpha, tau);
        let y_small = fractional_relaxation(y0, 1.0, alpha, tau);
        assert!(
            y_large < y_small,
            "y(10)={:.6} should be less than y(1)={:.6}",
            y_large,
            y_small
        );
    }

    #[test]
    fn test_fractional_oscillator_free_response() {
        let osc = FractionalOscillator::new(1.0, 0.0, 1.0, 0.8);
        // At t=0, x = x0
        let x0 = 2.0;
        let x = osc.free_response(x0, 0.0);
        assert!((x - x0).abs() < 1e-10);
        // Response should be bounded and start decaying
        let x1 = osc.free_response(x0, 1.0);
        assert!(x1.abs() <= x0.abs() + 1e-6);
    }

    #[test]
    fn test_fractional_oscillator_natural_freq() {
        let osc = FractionalOscillator::new(4.0, 0.0, 16.0, 1.0);
        // ω_n = sqrt(k/m) = sqrt(16/4) = 2
        assert!((osc.natural_frequency() - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_levy_stable_gaussian_limit() {
        // α=2 should produce Gaussian samples
        let params = LevyStableParams::gaussian(1.0);
        let samples: Vec<f64> = (0..1000)
            .map(|i| levy_stable_sample(&params, i as u64 * 7 + 13))
            .collect();
        let mean = samples.iter().sum::<f64>() / samples.len() as f64;
        let var = samples.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / samples.len() as f64;
        assert!(mean.abs() < 0.2, "mean = {:.6}", mean);
        assert!((var - 1.0).abs() < 0.5, "var = {:.6}", var);
    }

    #[test]
    fn test_levy_stable_pdf_normalizes() {
        // Integral of symmetric Lévy PDF ≈ 1 over a wide range
        let alpha = 1.5;
        let scale = 1.0;
        let n = 2000;
        let x_max = 30.0;
        let dx = 2.0 * x_max / n as f64;
        let mut integral = 0.0_f64;
        for k in 0..=n {
            let x = -x_max + k as f64 * dx;
            let p = levy_stable_pdf(x, alpha, scale, 500, 20.0);
            let w = if k == 0 || k == n { 0.5 } else { 1.0 };
            integral += w * p;
        }
        integral *= dx;
        assert!((integral - 1.0).abs() < 0.1, "integral = {:.6}", integral);
    }

    #[test]
    fn test_subdiffusion_trajectory_length() {
        let (x, y) = subdiffusion_trajectory(50, 0.5, 1.0, 42);
        assert_eq!(x.len(), 51);
        assert_eq!(y.len(), 51);
        assert_eq!(x[0], 0.0);
        assert_eq!(y[0], 0.0);
    }

    #[test]
    fn test_superdiffusion_trajectory_start() {
        let (x, y) = superdiffusion_trajectory(30, 1.5, 1.0, 99);
        assert_eq!(x.len(), 31);
        assert_eq!(y.len(), 31);
        assert_eq!(x[0], 0.0);
        assert_eq!(y[0], 0.0);
    }

    #[test]
    fn test_gl_derivative_integer_order() {
        // For α=1, GL derivative ~ first difference / h
        let h = 0.1;
        let n = 20;
        // f(t) = t  => f' = 1
        let f_vals: Vec<f64> = (0..=n).map(|k| k as f64 * h).collect();
        let d = rl_derivative(&f_vals, h, 1.0);
        // Should be approximately 1.0
        assert!((d - 1.0).abs() < 0.1, "GL d/dt[t] ≈ {:.6}", d);
    }

    #[test]
    fn test_fde_predictor_corrector_constant_rhs_matches_exact_solution() {
        // ^C D^0.5 y = 1, y(0) = 0  =>  y(t) = t^0.5 / Γ(1.5).  A constant
        // right-hand side is integrated exactly by the product-integration
        // weights of both the predictor and the corrector.
        let alpha = 0.5;
        let (t, y) = fde_predictor_corrector(|_t, _y| 1.0, 0.0, 1.0, 100, alpha);
        assert_eq!(t.len(), 101);
        assert_eq!(y.len(), 101);
        assert_eq!(y[0], 0.0);
        let expected = 1.0 / gamma(1.5);
        assert!(
            (y[100] - expected).abs() < 1e-8,
            "y(1) = {}, expected {expected}",
            y[100]
        );
        for (tk, yk) in t.iter().zip(y.iter()) {
            let exact = tk.powf(alpha) / gamma(alpha + 1.0);
            assert!((yk - exact).abs() < 1e-8, "t = {tk}: {yk} vs {exact}");
        }
    }

    #[test]
    fn test_fde_predictor_corrector_integer_order_is_exponential_decay() {
        // α = 1: y' = -y, y(0) = 1  =>  y(1) = e^{-1} (Heun's method).
        let (_t, y) = fde_predictor_corrector(|_t, y| -y, 1.0, 1.0, 100, 1.0);
        let err = (y[100] - (-1.0_f64).exp()).abs();
        assert!(err < 1e-4, "|y(1) - e^-1| = {err}");
        // Second-order convergence: halving h divides the error by about 4.
        let (_t, y2) = fde_predictor_corrector(|_t, y| -y, 1.0, 1.0, 200, 1.0);
        let err2 = (y2[200] - (-1.0_f64).exp()).abs();
        assert!(err / err2 > 3.5 && err / err2 < 4.5, "{err} / {err2}");
    }

    #[test]
    fn test_fde_predictor_corrector_constant_history_is_not_double_counted() {
        // With α = 1 and f ≡ 1 the solution is y0 + t.  The previous scheme
        // added y0 on top of the Grünwald-Letnikov memory and returned
        // y_1 = 2 y0 + h.
        let y0 = 3.0;
        let (t, y) = fde_predictor_corrector(|_t, _y| 1.0, y0, 2.0, 40, 1.0);
        for (tk, yk) in t.iter().zip(y.iter()) {
            assert!((yk - (y0 + tk)).abs() < 1e-12, "t = {tk}: y = {yk}");
        }
        // A zero right-hand side keeps the initial value for every order.
        for alpha in [0.3, 0.5, 0.9, 1.0] {
            let (_t, y) = fde_predictor_corrector(|_t, _y| 0.0, 2.0, 1.0, 25, alpha);
            assert!(
                y.iter().all(|&v| (v - 2.0).abs() < 1e-14),
                "alpha = {alpha}"
            );
        }
    }

    #[test]
    fn test_fde_predictor_corrector_fractional_relaxation() {
        // ^C D^0.5 y = -y, y(0) = 1  =>  y(t) = E_{0.5}(-t^0.5) = e^t erfc(√t),
        // so y(1) = e · erfc(1) = 0.42758357615580705.
        let reference = mittag_leffler_1(0.5, -1.0, 200);
        assert!((reference - 0.427_583_576_155_807_05).abs() < 1e-9);
        let errors: Vec<f64> = [100_usize, 200]
            .iter()
            .map(|&n| {
                let (_t, y) = fde_predictor_corrector(|_t, y| -y, 1.0, 1.0, n, 0.5);
                (y[n] - reference).abs()
            })
            .collect();
        assert!(errors[1] < 3e-5, "error at 200 steps: {}", errors[1]);
        // Order 1 + α = 1.5: halving h divides the error by 2^1.5 ≈ 2.8.
        assert!(errors[0] / errors[1] > 2.3, "{errors:?}");
    }

    #[test]
    fn test_bagley_torvik_no_force() {
        // With zero force and zero initial conditions, solution stays at zero
        let (_, x) = bagley_torvik(1.0, 0.1, 1.0, |_| 0.0, 0.0, 0.0, 1.0, 50);
        for &xi in &x {
            assert!(xi.abs() < 1e-12, "x={:.6e}", xi);
        }
    }

    #[test]
    fn test_fractional_laplacian_2d_shape() {
        let nx = 8;
        let ny = 8;
        let h = 0.1;
        let f: Vec<f64> = (0..nx * ny)
            .map(|i| {
                let ix = i % nx;
                let iy = i / nx;
                let x = (ix as f64 - nx as f64 / 2.0) * h;
                let y2 = (iy as f64 - ny as f64 / 2.0) * h;
                (-x * x - y2 * y2).exp()
            })
            .collect();
        let lap = fractional_laplacian_2d(&f, nx, ny, h, 0.5).expect("valid arguments");
        assert_eq!(lap.len(), nx * ny);
        for &v in &lap {
            assert!(v.is_finite());
        }
    }

    /// Gaussian bump `exp(-x² / 2)` sampled on the symmetric grid
    /// `x_i = (i - (n-1)/2) h` (so that `x_{n-1-i} = -x_i` exactly).
    fn symmetric_gaussian(n: usize, h: f64) -> Vec<f64> {
        (0..n)
            .map(|i| {
                let x = (i as f64 - (n as f64 - 1.0) / 2.0) * h;
                (-x * x / 2.0).exp()
            })
            .collect()
    }

    #[test]
    fn test_riesz_derivative_symmetric_input_gives_symmetric_output() {
        let n = 64;
        let f = symmetric_gaussian(n, 0.1);
        for alpha in [0.3, 0.5, 1.2, 1.5, 1.9] {
            let d = riesz_derivative(&f, 0.1, alpha);
            assert_eq!(d.len(), n);
            for j in 0..n {
                assert!(
                    (d[j] - d[n - 1 - j]).abs() < 1e-12,
                    "alpha = {alpha}, j = {j}: {} vs {}",
                    d[j],
                    d[n - 1 - j]
                );
            }
        }
    }

    /// Textbook unshifted Riesz-Grünwald-Letnikov sum written as plain index
    /// loops, sharing no code with [`riesz_derivative`] beyond the coefficients:
    /// `c_α h^{-α} [ Σ_{i≤j} g_{j-i} f_i + Σ_{i≥j} g_{i-j} f_i ]`.  The centre
    /// sample `f_j` enters both the left- and the right-sided sum.
    fn riesz_reference(f: &[f64], h: f64, alpha: f64) -> Vec<f64> {
        let n = f.len();
        let g = gl_coefficients(alpha, n);
        let c = -1.0 / (2.0 * (PI * alpha / 2.0).cos());
        let mut out = vec![0.0; n];
        for j in 0..n {
            let mut left = 0.0;
            for i in 0..=j {
                left += g[j - i] * f[i];
            }
            let mut right = 0.0;
            for i in j..n {
                right += g[i - j] * f[i];
            }
            out[j] = c * (left + right) / h.powf(alpha);
        }
        out
    }

    #[test]
    fn test_riesz_derivative_sums_run_over_the_whole_grid() {
        // Unit masses at both ends: the left-sided sum at x_j sees f_0 through
        // g_j and the right-sided sum sees f_{n-1} through g_{n-1-j}.  A
        // right-sided sum truncated at j + 1 terms misses f_{n-1} for j < n/2.
        let n = 12;
        let h = 0.5;
        for alpha in [0.5, 1.5] {
            let mut f = vec![0.0; n];
            f[0] = 1.0;
            f[n - 1] = 1.0;
            let d = riesz_derivative(&f, h, alpha);
            let expected = riesz_reference(&f, h, alpha);
            for (j, (dj, ej)) in d.iter().zip(expected.iter()).enumerate() {
                assert!(
                    (dj - ej).abs() < 1e-12 * ej.abs().max(1.0),
                    "alpha = {alpha}, j = {j}: {dj} vs {ej}"
                );
            }
            // The far-end mass alone must reach x_0 through g_{n-1}.
            let g = gl_coefficients(alpha, n);
            let c = -1.0 / (2.0 * (PI * alpha / 2.0).cos());
            let mut right_only = vec![0.0; n];
            right_only[n - 1] = 1.0;
            let d0 = riesz_derivative(&right_only, h, alpha)[0];
            assert!(
                (d0 - c * g[n - 1] / h.powf(alpha)).abs() < 1e-15,
                "alpha = {alpha}: d[0] = {d0}"
            );
            assert!(d0 != 0.0);
        }
    }

    #[test]
    fn test_riesz_derivative_matches_reference_sums_for_general_input() {
        // A non-symmetric signal exercises every term of both sums.
        let n = 30;
        let h = 0.2;
        let f: Vec<f64> = (0..n)
            .map(|i| (0.7 * i as f64).sin() + 0.1 * i as f64)
            .collect();
        for alpha in [0.3, 0.8, 1.2, 1.7] {
            let d = riesz_derivative(&f, h, alpha);
            let expected = riesz_reference(&f, h, alpha);
            assert_eq!(d.len(), n);
            for (j, (dj, ej)) in d.iter().zip(expected.iter()).enumerate() {
                assert!(
                    (dj - ej).abs() < 1e-11 * ej.abs().max(1.0),
                    "alpha = {alpha}, j = {j}: {dj} vs {ej}"
                );
            }
        }
    }

    #[test]
    fn test_riesz_derivative_converges_to_analytic_value() {
        // For f(x) = exp(-x²/2): ∂^α f/∂|x|^α at x = 0 equals
        // -2^{α/2} Γ((α+1)/2) / √π.  The GL scheme is first-order accurate.
        let alpha = 1.5;
        let exact = -(2.0_f64.powf(alpha / 2.0)) * gamma((alpha + 1.0) / 2.0) / PI.sqrt();
        let errors: Vec<f64> = [0.2, 0.1, 0.05]
            .iter()
            .map(|&h| {
                let n = (24.0_f64 / h).round() as usize + 1;
                let f: Vec<f64> = (0..n)
                    .map(|i| {
                        let x = -12.0 + i as f64 * h;
                        (-x * x / 2.0).exp()
                    })
                    .collect();
                let d = riesz_derivative(&f, h, alpha);
                (d[n / 2] - exact).abs() / exact.abs()
            })
            .collect();
        assert!(errors[2] < 0.07, "relative errors {errors:?}");
        for w in errors.windows(2) {
            let ratio = w[0] / w[1];
            assert!(ratio > 1.8 && ratio < 2.4, "first-order ratio {ratio}");
        }
    }

    #[test]
    fn test_riesz_derivative_empty_input() {
        assert!(riesz_derivative(&[], 0.1, 0.5).is_empty());
    }

    #[test]
    fn test_fractional_laplacian_1d_unit_order_is_second_difference() {
        let n = 17;
        let h = 0.3;
        let f: Vec<f64> = (0..n)
            .map(|i| (0.7 * i as f64).sin() + 0.1 * i as f64)
            .collect();
        let lap = fractional_laplacian_1d(&f, h, 1.0).expect("s = 1 is valid");
        assert_eq!(lap.len(), n);
        for j in 0..n {
            let left = if j > 0 { f[j - 1] } else { 0.0 };
            let right = if j + 1 < n { f[j + 1] } else { 0.0 };
            let expected = -(left - 2.0 * f[j] + right) / (h * h);
            assert!(
                (lap[j] - expected).abs() < 1e-12 * expected.abs().max(1.0),
                "j = {j}: {} vs {expected}",
                lap[j]
            );
        }
    }

    #[test]
    fn test_fractional_laplacian_1d_half_order_matches_lorentzian() {
        // f = 1/(1+x²) has Fourier transform π e^{-|ξ|}, so
        // (-Δ)^{1/2} f = F⁻¹[|ξ| π e^{-|ξ|}] = (1 - x²)/(1 + x²)².
        let h = 0.05;
        let n = 2001;
        let x: Vec<f64> = (0..n).map(|i| -50.0 + i as f64 * h).collect();
        let f: Vec<f64> = x.iter().map(|&xi| 1.0 / (1.0 + xi * xi)).collect();
        let lap = fractional_laplacian_1d(&f, h, 0.5).expect("s = 1/2 is valid");
        let mut checked = 0;
        for (xi, li) in x.iter().zip(lap.iter()) {
            if xi.abs() <= 5.0 {
                let exact = (1.0 - xi * xi) / (1.0 + xi * xi).powi(2);
                assert!((li - exact).abs() < 1.5e-3, "x = {xi}: {li} vs {exact}");
                checked += 1;
            }
        }
        assert!(checked > 100);
        // The operator is not identically zero at s = 1/2.
        assert!(lap.iter().any(|v| v.abs() > 0.1));
    }

    #[test]
    fn test_fractional_laplacian_1d_gaussian_centre_value() {
        // (-Δ)^s exp(-x²/2) at x = 0 equals 2^s Γ(s + 1/2) / √π.
        let h = 0.05;
        let n = 801;
        let f: Vec<f64> = (0..n)
            .map(|i| {
                let x = (i as f64 - 400.0) * h;
                (-x * x / 2.0).exp()
            })
            .collect();
        for s in [0.25, 0.5, 0.75, 1.0] {
            let lap = fractional_laplacian_1d(&f, h, s).expect("valid s");
            let exact = 2.0_f64.powf(s) * gamma(s + 0.5) / PI.sqrt();
            assert!(
                (lap[400] - exact).abs() < 2e-3,
                "s = {s}: {} vs {exact}",
                lap[400]
            );
        }
    }

    #[test]
    fn test_fractional_laplacian_1d_is_symmetric_positive_semidefinite() {
        let n = 40;
        let h = 0.25;
        let mut state = 0x1234_5678_u64;
        let mut next = move || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (state >> 11) as f64 / (1u64 << 53) as f64 - 0.5
        };
        let f: Vec<f64> = (0..n).map(|_| next()).collect();
        let g: Vec<f64> = (0..n).map(|_| next()).collect();
        for s in [0.3, 0.5, 0.8, 1.0] {
            let lf = fractional_laplacian_1d(&f, h, s).expect("valid s");
            let lg = fractional_laplacian_1d(&g, h, s).expect("valid s");
            let quad: f64 = f.iter().zip(lf.iter()).map(|(a, b)| a * b).sum();
            assert!(quad > 0.0, "s = {s}: <f, Lf> = {quad}");
            let fg: f64 = f.iter().zip(lg.iter()).map(|(a, b)| a * b).sum();
            let gf: f64 = g.iter().zip(lf.iter()).map(|(a, b)| a * b).sum();
            assert!((fg - gf).abs() < 1e-12 * fg.abs().max(1.0), "symmetry");
        }
    }

    #[test]
    fn test_fractional_laplacian_1d_single_point_and_empty() {
        let lap = fractional_laplacian_1d(&[2.0], 0.5, 0.5).expect("valid");
        let expected = 2.0 * gamma(2.0) / (gamma(1.5) * gamma(1.5)) / 0.5;
        assert!(
            (lap[0] - expected).abs() < 1e-12,
            "{} vs {expected}",
            lap[0]
        );
        assert!(
            fractional_laplacian_1d(&[], 0.5, 0.5)
                .expect("valid")
                .is_empty()
        );
    }

    #[test]
    fn test_fractional_laplacian_rejects_invalid_arguments() {
        let f = vec![1.0, 2.0, 3.0, 4.0];
        for s in [0.0, -0.5, 1.5, f64::NAN, f64::INFINITY] {
            match fractional_laplacian_1d(&f, 0.1, s) {
                Err(PhysicsError::InvalidInput { field, .. }) => assert_eq!(field, "s"),
                other => panic!("s = {s}: expected InvalidInput, got {other:?}"),
            }
        }
        for h in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            match fractional_laplacian_1d(&f, h, 0.5) {
                Err(PhysicsError::InvalidInput { field, .. }) => assert_eq!(field, "h"),
                other => panic!("h = {h}: expected InvalidInput, got {other:?}"),
            }
        }
        // 2-D: invalid s, invalid h, and a data length that does not match the grid.
        assert!(fractional_laplacian_2d(&f, 2, 2, 0.1, 1.2).is_err());
        assert!(fractional_laplacian_2d(&f, 2, 2, 0.0, 0.5).is_err());
        match fractional_laplacian_2d(&f, 3, 2, 0.1, 0.5) {
            Err(PhysicsError::InvalidInput { field, .. }) => assert_eq!(field, "f"),
            other => panic!("expected InvalidInput, got {other:?}"),
        }
        assert!(fractional_laplacian_2d(&f, usize::MAX, 2, 0.1, 0.5).is_err());
    }

    #[test]
    fn test_fractional_laplacian_2d_unit_order_is_five_point_laplacian() {
        let nx = 7;
        let ny = 5;
        let h = 0.4;
        let f: Vec<f64> = (0..nx * ny)
            .map(|i| ((i * 37 % 11) as f64 - 5.0) * 0.3 + (i as f64 * 0.21).cos())
            .collect();
        let lap = fractional_laplacian_2d(&f, nx, ny, h, 1.0).expect("s = 1 is valid");
        let at = |ix: isize, iy: isize| -> f64 {
            if ix < 0 || iy < 0 || ix >= nx as isize || iy >= ny as isize {
                0.0
            } else {
                f[iy as usize * nx + ix as usize]
            }
        };
        for iy in 0..ny as isize {
            for ix in 0..nx as isize {
                let expected = -(at(ix - 1, iy) + at(ix + 1, iy) + at(ix, iy - 1) + at(ix, iy + 1)
                    - 4.0 * at(ix, iy))
                    / (h * h);
                let got = lap[iy as usize * nx + ix as usize];
                assert!(
                    (got - expected).abs() < 1e-11 * expected.abs().max(1.0),
                    "({ix}, {iy}): {got} vs {expected}"
                );
            }
        }
    }

    #[test]
    fn test_fractional_laplacian_2d_is_sum_of_axis_operators() {
        // For a separable field u(x) v(y) the axis-wise operator is
        // (L u)(x) v(y) + u(x) (L v)(y) with L the 1-D fractional Laplacian.
        let nx = 9;
        let ny = 6;
        let h = 0.3;
        let s = 0.5;
        let u: Vec<f64> = (0..nx).map(|i| (0.5 * i as f64).sin() + 1.0).collect();
        let v: Vec<f64> = (0..ny).map(|j| (0.8 * j as f64).cos() + 0.5).collect();
        let f: Vec<f64> = (0..nx * ny).map(|k| u[k % nx] * v[k / nx]).collect();
        let lap = fractional_laplacian_2d(&f, nx, ny, h, s).expect("valid");
        let lu = fractional_laplacian_1d(&u, h, s).expect("valid");
        let lv = fractional_laplacian_1d(&v, h, s).expect("valid");
        for (k, got) in lap.iter().enumerate() {
            let (ix, iy) = (k % nx, k / nx);
            let expected = lu[ix] * v[iy] + u[ix] * lv[iy];
            assert!(
                (got - expected).abs() < 1e-12 * expected.abs().max(1.0),
                "k = {k}: {got} vs {expected}"
            );
        }
        // Degenerate grids.
        assert!(
            fractional_laplacian_2d(&[], 0, 5, h, s)
                .expect("valid")
                .is_empty()
        );
        assert!(
            fractional_laplacian_2d(&[], 4, 0, h, s)
                .expect("valid")
                .is_empty()
        );
    }

    #[test]
    fn test_memory_kernel_mittag_leffler_at_zero() {
        // K(0, alpha, lambda) = 0 by definition (t <= 0)
        assert_eq!(memory_kernel_mittag_leffler(0.0, 0.5, 1.0), 0.0);
    }

    #[test]
    fn test_relaxation_time_threshold() {
        // Exponential relaxation α=1 reaches threshold faster than α=0.3
        let t1 = relaxation_time_threshold(1.0, 1.0, 0.01, 20.0, 1000);
        let t2 = relaxation_time_threshold(0.3, 1.0, 0.01, 200.0, 10000);
        assert!(t1.is_some());
        assert!(t2.is_some());
        assert!(
            t1.unwrap() < t2.unwrap(),
            "α=1 reaches threshold at t={:.4}, α=0.3 at t={:.4}",
            t1.unwrap(),
            t2.unwrap()
        );
    }

    #[test]
    fn test_log_gamma_positive() {
        // ln Γ(5) = ln 24
        let lg = log_gamma(5.0);
        assert!((lg - 24.0_f64.ln()).abs() < 1e-10);
    }

    #[test]
    fn test_mittag_leffler_derivative() {
        // d/dz E_α(z) = E_{α,α}(z)
        let alpha = 0.7;
        let z = -0.5;
        let deriv = mittag_leffler_derivative(alpha, z, 100);
        let direct = mittag_leffler_2(alpha, alpha, z, 100);
        assert!((deriv - direct).abs() < 1e-12);
    }

    #[test]
    fn test_rl_integral_series_length() {
        let f = vec![1.0; 10];
        let series = rl_integral_series(&f, 0.1, 0.5);
        assert_eq!(series.len(), 10);
    }

    #[test]
    fn test_rl_derivative_series_length() {
        let f = vec![1.0; 10];
        let series = rl_derivative_series(&f, 0.1, 0.5);
        assert_eq!(series.len(), 10);
    }

    #[test]
    fn test_caputo_series_length() {
        let f: Vec<f64> = (0..20).map(|k| k as f64 * 0.1).collect();
        let series = caputo_derivative_series(&f, 0.1, 0.5);
        assert_eq!(series.len(), 19);
    }

    #[test]
    fn test_subdiffusion_pdf_normalizes_approx() {
        // PDF should integrate to ≈ 1 over a wide range
        let params = AnomalousDiffusionParams::new(0.5, 0.1);
        let t = 1.0;
        let n = 1000;
        let x_max = 5.0;
        let dx = 2.0 * x_max / n as f64;
        let integral: f64 = (0..=n)
            .map(|k| {
                let x = -x_max + k as f64 * dx;
                let p = subdiffusion_pdf(x, t, &params);
                let w = if k == 0 || k == n { 0.5 } else { 1.0 };
                w * p
            })
            .sum::<f64>()
            * dx;
        assert!((integral - 1.0).abs() < 0.05, "integral = {:.6}", integral);
    }
}
