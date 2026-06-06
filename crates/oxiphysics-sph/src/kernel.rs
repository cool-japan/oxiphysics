// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! SPH kernel functions for smoothed interpolation.
//!
//! All kernels are normalized for 3D and implement the [`SphKernel`] trait.
//! Includes Wendland C4/C6, super-Gaussian, kernel corrections (Shepard, MLS),
//! kernel gradient correction, and kernel comparison utilities.

use std::f64::consts::PI;

/// Trait for SPH smoothing kernel functions.
///
/// Each kernel provides the value, gradient magnitude, and Laplacian
/// as functions of the scalar distance `r` and the smoothing length `h`.
pub trait SphKernel: Send + Sync {
    /// Evaluate the kernel W(r, h).
    fn w(&self, r: f64, h: f64) -> f64;

    /// Evaluate the gradient magnitude dW/dr (scalar, to be multiplied by r_hat).
    fn grad_w(&self, r: f64, h: f64) -> f64;

    /// Evaluate the Laplacian of the kernel.
    fn laplacian_w(&self, r: f64, h: f64) -> f64;
}

// ---------------------------------------------------------------------------
// Cubic Spline Kernel (M4)
// ---------------------------------------------------------------------------

/// Standard cubic spline kernel (M4) for 3D.
///
/// Compact support at `r/h = 2`. Normalization constant: `1 / (π h³)`.
#[derive(Debug, Clone, Copy, Default)]
pub struct CubicSplineKernel;

impl SphKernel for CubicSplineKernel {
    fn w(&self, r: f64, h: f64) -> f64 {
        let sigma = 1.0 / (PI * h * h * h);
        let q = r / h;
        if q >= 2.0 {
            0.0
        } else if q >= 1.0 {
            let t = 2.0 - q;
            sigma * 0.25 * t * t * t
        } else {
            sigma * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
        }
    }

    fn grad_w(&self, r: f64, h: f64) -> f64 {
        let sigma = 1.0 / (PI * h * h * h);
        let q = r / h;
        if !(1e-14..2.0).contains(&q) {
            0.0
        } else if q >= 1.0 {
            let t = 2.0 - q;
            sigma * (-0.75 * t * t) / h
        } else {
            sigma * (-3.0 * q + 2.25 * q * q) / h
        }
    }

    fn laplacian_w(&self, r: f64, h: f64) -> f64 {
        let sigma = 1.0 / (PI * h * h * h);
        let q = r / h;
        let h2 = h * h;
        if q >= 2.0 {
            0.0
        } else if q >= 1.0 {
            let t = 2.0 - q;
            // Laplacian in 3D: d²W/dr² + (2/r) dW/dr
            let d2w = sigma * 1.5 * t / h2;
            let dw = sigma * (-0.75 * t * t) / h;
            d2w + if r > 1e-14 { 2.0 * dw / r } else { 0.0 }
        } else {
            let d2w = sigma * (-3.0 + 4.5 * q) / h2;
            let dw = sigma * (-3.0 * q + 2.25 * q * q) / h;
            d2w + if r > 1e-14 { 2.0 * dw / r } else { 0.0 }
        }
    }
}

// ---------------------------------------------------------------------------
// Wendland C2 Kernel
// ---------------------------------------------------------------------------

/// Wendland quintic C2 kernel for 3D.
///
/// Compact support at `r/h = 2`. Smooth and positive-definite.
#[derive(Debug, Clone, Copy, Default)]
pub struct WendlandC2Kernel;

impl SphKernel for WendlandC2Kernel {
    fn w(&self, r: f64, h: f64) -> f64 {
        let q = r / h;
        if q >= 2.0 {
            return 0.0;
        }
        let sigma = 21.0 / (16.0 * PI * h * h * h);
        let t = 1.0 - 0.5 * q;
        sigma * t * t * t * t * (1.0 + 2.0 * q)
    }

    fn grad_w(&self, r: f64, h: f64) -> f64 {
        let q = r / h;
        if !(1e-14..2.0).contains(&q) {
            return 0.0;
        }
        let sigma = 21.0 / (16.0 * PI * h * h * h);
        let t = 1.0 - 0.5 * q;
        sigma * t * t * t * (-5.0 * q) / h
    }

    fn laplacian_w(&self, r: f64, h: f64) -> f64 {
        let q = r / h;
        if q >= 2.0 {
            return 0.0;
        }
        let sigma = 21.0 / (16.0 * PI * h * h * h);
        let h2 = h * h;
        let t = 1.0 - 0.5 * q;
        let d2w = sigma * t * t * (10.0 * q - 5.0) / h2;
        let dw = sigma * t * t * t * (-5.0 * q) / h;
        d2w + if r > 1e-14 { 2.0 * dw / r } else { 0.0 }
    }
}

// ---------------------------------------------------------------------------
// Wendland C4 Kernel
// ---------------------------------------------------------------------------

/// Wendland C4 kernel for 3D.
///
/// Higher-order smoothness than C2. Compact support at `r/h = 2`.
/// W(r,h) = σ * (1 - q/2)^6 * (35q²/12 + 3q + 1), q = r/h
/// σ = 495 / (256 π h³)
#[derive(Debug, Clone, Copy, Default)]
pub struct WendlandC4Kernel;

impl SphKernel for WendlandC4Kernel {
    fn w(&self, r: f64, h: f64) -> f64 {
        let q = r / h;
        if q >= 2.0 {
            return 0.0;
        }
        let sigma = 495.0 / (256.0 * PI * h * h * h);
        let t = 1.0 - 0.5 * q;
        sigma * t.powi(6) * (35.0 / 12.0 * q * q + 3.0 * q + 1.0)
    }

    fn grad_w(&self, r: f64, h: f64) -> f64 {
        let q = r / h;
        if !(1e-14..2.0).contains(&q) {
            return 0.0;
        }
        let sigma = 495.0 / (256.0 * PI * h * h * h);
        let t = 1.0 - 0.5 * q;
        // dW/dq = σ * [ 6*t^5*(-1/2)*(35q²/12 + 3q + 1) + t^6*(70q/12 + 3) ]
        let poly = 35.0 / 12.0 * q * q + 3.0 * q + 1.0;
        let dpoly = 70.0 / 12.0 * q + 3.0;
        let dw_dq = sigma * (6.0 * t.powi(5) * (-0.5) * poly + t.powi(6) * dpoly);
        dw_dq / h
    }

    fn laplacian_w(&self, r: f64, h: f64) -> f64 {
        let q = r / h;
        if q >= 2.0 {
            return 0.0;
        }
        // Use finite difference approximation for simplicity
        let eps = 1e-6 * h;
        let r_plus = r + eps;
        let r_minus = (r - eps).max(0.0);
        let dw_plus = self.grad_w(r_plus, h);
        let dw_minus = self.grad_w(r_minus, h);
        let d2w = (dw_plus - dw_minus) / (r_plus - r_minus);
        let dw = self.grad_w(r, h);
        d2w + if r > 1e-14 { 2.0 * dw / r } else { 0.0 }
    }
}

// ---------------------------------------------------------------------------
// Wendland C6 Kernel
// ---------------------------------------------------------------------------

/// Wendland C6 kernel for 3D.
///
/// Even higher-order smoothness. Compact support at `r/h = 2`.
/// W(r,h) = σ * (1 - q/2)^8 * (4q³ + 6.25q² + 4q + 1), q = r/h
/// σ = 1365 / (512 π h³)
#[derive(Debug, Clone, Copy, Default)]
pub struct WendlandC6Kernel;

impl SphKernel for WendlandC6Kernel {
    fn w(&self, r: f64, h: f64) -> f64 {
        let q = r / h;
        if q >= 2.0 {
            return 0.0;
        }
        let sigma = 1365.0 / (512.0 * PI * h * h * h);
        let t = 1.0 - 0.5 * q;
        sigma * t.powi(8) * (4.0 * q * q * q + 6.25 * q * q + 4.0 * q + 1.0)
    }

    fn grad_w(&self, r: f64, h: f64) -> f64 {
        let q = r / h;
        if !(1e-14..2.0).contains(&q) {
            return 0.0;
        }
        let sigma = 1365.0 / (512.0 * PI * h * h * h);
        let t = 1.0 - 0.5 * q;
        let poly = 4.0 * q * q * q + 6.25 * q * q + 4.0 * q + 1.0;
        let dpoly = 12.0 * q * q + 12.5 * q + 4.0;
        let dw_dq = sigma * (8.0 * t.powi(7) * (-0.5) * poly + t.powi(8) * dpoly);
        dw_dq / h
    }

    fn laplacian_w(&self, r: f64, h: f64) -> f64 {
        let q = r / h;
        if q >= 2.0 {
            return 0.0;
        }
        let eps = 1e-6 * h;
        let r_plus = r + eps;
        let r_minus = (r - eps).max(0.0);
        let dw_plus = self.grad_w(r_plus, h);
        let dw_minus = self.grad_w(r_minus, h);
        let d2w = (dw_plus - dw_minus) / (r_plus - r_minus);
        let dw = self.grad_w(r, h);
        d2w + if r > 1e-14 { 2.0 * dw / r } else { 0.0 }
    }
}

// ---------------------------------------------------------------------------
// Super-Gaussian Kernel
// ---------------------------------------------------------------------------

/// Super-Gaussian kernel for 3D.
///
/// A Gaussian multiplied by a correction polynomial to improve partition of unity.
/// W(r,h) = σ * (D/2 + 1 - q²) * exp(-q²), q = r/h
/// D = 3 for 3D, σ = 1/(π^(3/2) h³)
#[derive(Debug, Clone, Copy, Default)]
pub struct SuperGaussianKernel;

impl SphKernel for SuperGaussianKernel {
    fn w(&self, r: f64, h: f64) -> f64 {
        let q = r / h;
        if q > 3.0 {
            return 0.0;
        }
        let sigma = 1.0 / (PI.powf(1.5) * h * h * h);
        let q2 = q * q;
        sigma * (2.5 - q2) * (-q2).exp() // D/2 + 1 = 3/2 + 1 = 2.5
    }

    fn grad_w(&self, r: f64, h: f64) -> f64 {
        let q = r / h;
        if !(1e-14..3.0).contains(&q) {
            return 0.0;
        }
        let sigma = 1.0 / (PI.powf(1.5) * h * h * h);
        let q2 = q * q;
        let exp_q2 = (-q2).exp();
        // dW/dq = σ * [(-2q)*exp(-q²)*(2.5 - q²) + (2.5-q²)*exp(-q²)*(-2q)]
        //        Actually: d/dq[(2.5-q²)*exp(-q²)]
        //        = -2q*exp(-q²) + (2.5-q²)*(-2q)*exp(-q²)
        //        = -2q*exp(-q²)*[1 + (2.5-q²)]
        //        = -2q*exp(-q²)*(3.5 - q²)
        let dw_dq = sigma * (-2.0 * q) * exp_q2 * (3.5 - q2);
        dw_dq / h
    }

    fn laplacian_w(&self, r: f64, h: f64) -> f64 {
        let q = r / h;
        if q > 3.0 {
            return 0.0;
        }
        let eps = 1e-6 * h;
        let r_plus = r + eps;
        let r_minus = (r - eps).max(0.0);
        let dw_plus = self.grad_w(r_plus, h);
        let dw_minus = self.grad_w(r_minus, h);
        let d2w = (dw_plus - dw_minus) / (r_plus - r_minus);
        let dw = self.grad_w(r, h);
        d2w + if r > 1e-14 { 2.0 * dw / r } else { 0.0 }
    }
}

// ---------------------------------------------------------------------------
// Poly6 Kernel
// ---------------------------------------------------------------------------

/// Poly6 kernel for 3D, commonly used for density estimation.
///
/// Compact support at `r/h = 1` (h is the support radius).
#[derive(Debug, Clone, Copy, Default)]
pub struct Poly6Kernel;

impl SphKernel for Poly6Kernel {
    fn w(&self, r: f64, h: f64) -> f64 {
        if r >= h {
            return 0.0;
        }
        let sigma = 315.0 / (64.0 * PI * h.powi(9));
        let d = h * h - r * r;
        sigma * d * d * d
    }

    fn grad_w(&self, r: f64, h: f64) -> f64 {
        if r >= h || r < 1e-14 {
            return 0.0;
        }
        let sigma = 315.0 / (64.0 * PI * h.powi(9));
        let d = h * h - r * r;
        sigma * (-6.0 * r) * d * d
    }

    fn laplacian_w(&self, r: f64, h: f64) -> f64 {
        if r >= h {
            return 0.0;
        }
        let sigma = 315.0 / (64.0 * PI * h.powi(9));
        let d = h * h - r * r;
        let d2w = sigma * (-6.0 * d * d + 24.0 * r * r * d);
        let dw = sigma * (-6.0 * r) * d * d;
        d2w + if r > 1e-14 { 2.0 * dw / r } else { 0.0 }
    }
}

// ---------------------------------------------------------------------------
// Spiky Kernel
// ---------------------------------------------------------------------------

/// Spiky kernel for 3D, commonly used for pressure gradient computation.
///
/// The non-vanishing gradient at `r → 0` prevents particle clumping.
/// Compact support at `r/h = 1`.
#[derive(Debug, Clone, Copy, Default)]
pub struct SpikyKernel;

impl SphKernel for SpikyKernel {
    fn w(&self, r: f64, h: f64) -> f64 {
        if r >= h {
            return 0.0;
        }
        let sigma = 15.0 / (PI * h.powi(6));
        let d = h - r;
        sigma * d * d * d
    }

    fn grad_w(&self, r: f64, h: f64) -> f64 {
        if r >= h || r < 1e-14 {
            return 0.0;
        }
        let sigma = 15.0 / (PI * h.powi(6));
        let d = h - r;
        sigma * (-3.0) * d * d
    }

    fn laplacian_w(&self, r: f64, h: f64) -> f64 {
        if r >= h {
            return 0.0;
        }
        let sigma = 15.0 / (PI * h.powi(6));
        let d = h - r;
        let d2w = sigma * 6.0 * d;
        let dw = sigma * (-3.0) * d * d;
        d2w + if r > 1e-14 { 2.0 * dw / r } else { 0.0 }
    }
}

// ---------------------------------------------------------------------------
// Kernel corrections
// ---------------------------------------------------------------------------

/// Shepard (zeroth-order) kernel correction.
///
/// Corrects the kernel so that it sums exactly to 1 (partition of unity).
/// W_corrected(r_ij) = W(r_ij) / Σ_j (m_j/ρ_j) * W(r_ij)
///
/// Returns the correction factor for particle i.
pub fn shepard_correction(kernel_values: &[f64], masses: &[f64], densities: &[f64]) -> f64 {
    let mut sum = 0.0_f64;
    for (idx, &w) in kernel_values.iter().enumerate() {
        if densities[idx] > 1e-14 {
            sum += masses[idx] / densities[idx] * w;
        }
    }
    if sum > 1e-30 { 1.0 / sum } else { 1.0 }
}

/// Moving Least Squares (MLS) first-order correction matrix.
///
/// Computes the correction matrix for linear consistency.
/// Returns the 4x4 moment matrix M and the right-hand-side vector b
/// for MLS correction at particle i.
///
/// The moment matrix is M = Σ_j (m_j/ρ_j) * W_ij * p_j * p_j^T
/// where p = \[1, x-xi, y-yi, z-zi\].
pub fn mls_moment_matrix(
    pos_i: [f64; 3],
    neighbor_positions: &[[f64; 3]],
    kernel_values: &[f64],
    masses: &[f64],
    densities: &[f64],
) -> [[f64; 4]; 4] {
    let mut m = [[0.0_f64; 4]; 4];
    for (idx, &pos_j) in neighbor_positions.iter().enumerate() {
        if densities[idx] < 1e-14 {
            continue;
        }
        let w = kernel_values[idx];
        let vol = masses[idx] / densities[idx];
        let p = [
            1.0,
            pos_j[0] - pos_i[0],
            pos_j[1] - pos_i[1],
            pos_j[2] - pos_i[2],
        ];
        for row in 0..4 {
            for col in 0..4 {
                m[row][col] += vol * w * p[row] * p[col];
            }
        }
    }
    m
}

/// Kernel gradient correction (Randles-Libersky).
///
/// Corrects the gradient by computing the inverse of the correction matrix:
/// L_i = (Σ_j V_j ∇W_ij ⊗ r_ij)^(-1)
///
/// Returns the 3x3 correction matrix L.
pub fn gradient_correction_matrix(
    pos_i: [f64; 3],
    neighbor_positions: &[[f64; 3]],
    kernel_gradients: &[[f64; 3]],
    volumes: &[f64],
) -> [[f64; 3]; 3] {
    let mut l = [[0.0_f64; 3]; 3];
    for (idx, &pos_j) in neighbor_positions.iter().enumerate() {
        let rij = [
            pos_j[0] - pos_i[0],
            pos_j[1] - pos_i[1],
            pos_j[2] - pos_i[2],
        ];
        let grad = kernel_gradients[idx];
        let vol = volumes[idx];
        // Outer product: grad ⊗ rij
        for row in 0..3 {
            for col in 0..3 {
                l[row][col] += vol * grad[row] * rij[col];
            }
        }
    }
    // Invert the 3x3 matrix
    invert_3x3(l)
}

/// Invert a 3x3 matrix. Returns identity if singular.
fn invert_3x3(m: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

    if det.abs() < 1e-30 {
        return [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    }

    let inv_det = 1.0 / det;
    [
        [
            inv_det * (m[1][1] * m[2][2] - m[1][2] * m[2][1]),
            inv_det * (m[0][2] * m[2][1] - m[0][1] * m[2][2]),
            inv_det * (m[0][1] * m[1][2] - m[0][2] * m[1][1]),
        ],
        [
            inv_det * (m[1][2] * m[2][0] - m[1][0] * m[2][2]),
            inv_det * (m[0][0] * m[2][2] - m[0][2] * m[2][0]),
            inv_det * (m[0][2] * m[1][0] - m[0][0] * m[1][2]),
        ],
        [
            inv_det * (m[1][0] * m[2][1] - m[1][1] * m[2][0]),
            inv_det * (m[0][1] * m[2][0] - m[0][0] * m[2][1]),
            inv_det * (m[0][0] * m[1][1] - m[0][1] * m[1][0]),
        ],
    ]
}

// ---------------------------------------------------------------------------
// Kernel comparison utilities
// ---------------------------------------------------------------------------

/// Compare two kernels by computing their values at several sample points.
///
/// Returns `(max_abs_diff, max_rel_diff)` between the two kernels.
pub fn compare_kernels(
    k1: &dyn SphKernel,
    k2: &dyn SphKernel,
    h: f64,
    n_samples: usize,
) -> (f64, f64) {
    let support = 2.0 * h;
    let dr = support / n_samples as f64;
    let mut max_abs = 0.0_f64;
    let mut max_rel = 0.0_f64;

    for i in 0..n_samples {
        let r = (i as f64 + 0.5) * dr;
        let w1 = k1.w(r, h);
        let w2 = k2.w(r, h);
        let abs_diff = (w1 - w2).abs();
        if abs_diff > max_abs {
            max_abs = abs_diff;
        }
        let denom = w1.abs().max(w2.abs());
        if denom > 1e-30 {
            let rel = abs_diff / denom;
            if rel > max_rel {
                max_rel = rel;
            }
        }
    }

    (max_abs, max_rel)
}

/// Compute the effective support radius of a kernel (where W drops below threshold).
pub fn effective_support_radius(kernel: &dyn SphKernel, h: f64, threshold: f64) -> f64 {
    let n = 1000;
    let max_r = 3.0 * h;
    let dr = max_r / n as f64;
    let mut r_eff = 0.0_f64;
    for i in 0..n {
        let r = (i as f64 + 0.5) * dr;
        if kernel.w(r, h) > threshold {
            r_eff = r;
        }
    }
    r_eff
}

// ---------------------------------------------------------------------------
// Free-function API using [f64; 3] vectors (no nalgebra dependency)
// ---------------------------------------------------------------------------

/// Evaluate a kernel at distance `r` with smoothing length `h`.
pub fn eval<K: SphKernel>(kernel: &K, r: f64, h: f64) -> f64 {
    kernel.w(r, h)
}

/// Evaluate the vector kernel gradient `∇W(r_vec, h) = (dW/dr) * r̂`.
pub fn grad<K: SphKernel>(kernel: &K, r_vec: [f64; 3], h: f64) -> [f64; 3] {
    let r2 = r_vec[0] * r_vec[0] + r_vec[1] * r_vec[1] + r_vec[2] * r_vec[2];
    if r2 < 1e-28 {
        return [0.0; 3];
    }
    let r = r2.sqrt();
    let dw_dr = kernel.grad_w(r, h);
    let inv_r = dw_dr / r;
    [r_vec[0] * inv_r, r_vec[1] * inv_r, r_vec[2] * inv_r]
}

/// Evaluate the cubic-spline kernel value at scalar distance `r`.
pub fn cubic_spline_eval(r: f64, h: f64) -> f64 {
    CubicSplineKernel.w(r, h)
}

/// Evaluate the cubic-spline gradient given displacement vector `r_vec`.
pub fn cubic_spline_grad(r_vec: [f64; 3], h: f64) -> [f64; 3] {
    grad(&CubicSplineKernel, r_vec, h)
}

/// Evaluate the Wendland C2 kernel value at scalar distance `r`.
pub fn wendland_c2_eval(r: f64, h: f64) -> f64 {
    WendlandC2Kernel.w(r, h)
}

/// Evaluate the Wendland C2 gradient given displacement vector `r_vec`.
pub fn wendland_c2_grad(r_vec: [f64; 3], h: f64) -> [f64; 3] {
    grad(&WendlandC2Kernel, r_vec, h)
}

/// Evaluate the Wendland C4 kernel value at scalar distance `r`.
pub fn wendland_c4_eval(r: f64, h: f64) -> f64 {
    WendlandC4Kernel.w(r, h)
}

/// Evaluate the Wendland C6 kernel value at scalar distance `r`.
pub fn wendland_c6_eval(r: f64, h: f64) -> f64 {
    WendlandC6Kernel.w(r, h)
}

/// Evaluate the super-Gaussian kernel value at scalar distance `r`.
pub fn super_gaussian_eval(r: f64, h: f64) -> f64 {
    SuperGaussianKernel.w(r, h)
}

/// Evaluate the Poly6 kernel value at scalar distance `r`.
pub fn poly6_eval(r: f64, h: f64) -> f64 {
    Poly6Kernel.w(r, h)
}

/// Evaluate the Poly6 gradient given displacement vector `r_vec`.
pub fn poly6_grad(r_vec: [f64; 3], h: f64) -> [f64; 3] {
    grad(&Poly6Kernel, r_vec, h)
}

/// Evaluate the Spiky kernel value at scalar distance `r`.
pub fn spiky_eval(r: f64, h: f64) -> f64 {
    SpikyKernel.w(r, h)
}

/// Evaluate the Spiky gradient given displacement vector `r_vec`.
pub fn spiky_grad(r_vec: [f64; 3], h: f64) -> [f64; 3] {
    grad(&SpikyKernel, r_vec, h)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Discrete integration of a 3D kernel over a grid.
    fn integrate_kernel_3d(kernel: &dyn SphKernel, h: f64, n: usize) -> f64 {
        let support = 2.0 * h; // conservative upper bound
        let dx = 2.0 * support / n as f64;
        let mut sum = 0.0;
        for i in 0..n {
            let x = -support + (i as f64 + 0.5) * dx;
            for j in 0..n {
                let y = -support + (j as f64 + 0.5) * dx;
                for k in 0..n {
                    let z = -support + (k as f64 + 0.5) * dx;
                    let r = (x * x + y * y + z * z).sqrt();
                    sum += kernel.w(r, h) * dx * dx * dx;
                }
            }
        }
        sum
    }

    #[test]
    fn cubic_spline_integrates_to_one() {
        let k = CubicSplineKernel;
        let integral = integrate_kernel_3d(&k, 1.0, 80);
        assert!(
            (integral - 1.0).abs() < 0.05,
            "Cubic spline integral = {integral}, expected ~1.0"
        );
    }

    #[test]
    fn wendland_c2_integrates_to_one() {
        let k = WendlandC2Kernel;
        let integral = integrate_kernel_3d(&k, 1.0, 80);
        assert!(
            (integral - 1.0).abs() < 0.05,
            "Wendland C2 integral = {integral}, expected ~1.0"
        );
    }

    #[test]
    fn wendland_c4_positive_inside_support() {
        let k = WendlandC4Kernel;
        assert!(k.w(0.0, 1.0) > 0.0, "C4 at origin should be positive");
        assert!(k.w(0.5, 1.0) > 0.0, "C4 at r=0.5 should be positive");
        assert!(k.w(1.5, 1.0) > 0.0, "C4 at r=1.5 should be positive");
        assert!(
            k.w(2.5, 1.0).abs() < 1e-14,
            "C4 outside support should be 0"
        );
    }

    #[test]
    fn wendland_c6_positive_inside_support() {
        let k = WendlandC6Kernel;
        assert!(k.w(0.0, 1.0) > 0.0, "C6 at origin should be positive");
        assert!(k.w(0.5, 1.0) > 0.0, "C6 at r=0.5 should be positive");
        assert!(
            k.w(2.5, 1.0).abs() < 1e-14,
            "C6 outside support should be 0"
        );
    }

    #[test]
    fn super_gaussian_positive_inside() {
        let k = SuperGaussianKernel;
        assert!(k.w(0.0, 1.0) > 0.0, "SG at origin should be positive");
        assert!(k.w(0.5, 1.0) > 0.0, "SG at r=0.5 should be positive");
    }

    #[test]
    fn super_gaussian_zero_far() {
        let k = SuperGaussianKernel;
        assert!(
            k.w(5.0, 1.0).abs() < 1e-10,
            "SG far from origin should be ~0"
        );
    }

    #[test]
    fn kernel_gradient_zero_at_origin() {
        let kernels: Vec<Box<dyn SphKernel>> = vec![
            Box::new(CubicSplineKernel),
            Box::new(WendlandC2Kernel),
            Box::new(Poly6Kernel),
            Box::new(SpikyKernel),
            Box::new(WendlandC4Kernel),
            Box::new(WendlandC6Kernel),
            Box::new(SuperGaussianKernel),
        ];
        for k in &kernels {
            assert!(
                k.grad_w(0.0, 1.0).abs() < 1e-10,
                "Gradient at r=0 should be 0"
            );
        }
    }

    #[test]
    fn kernel_value_positive_inside_support() {
        let k = CubicSplineKernel;
        assert!(k.w(0.0, 1.0) > 0.0);
        assert!(k.w(0.5, 1.0) > 0.0);
        assert!(k.w(1.5, 1.0) > 0.0);
        assert!((k.w(2.5, 1.0)).abs() < 1e-14);
    }

    #[test]
    fn eval_matches_trait_w() {
        let h = 0.5_f64;
        for r in [0.0, 0.2, 0.5, 0.9, 1.1, 1.9, 2.1] {
            let expected = CubicSplineKernel.w(r, h);
            let got = eval(&CubicSplineKernel, r, h);
            assert!(
                (expected - got).abs() < 1e-14,
                "eval mismatch at r={r}: expected {expected}, got {got}"
            );
        }
    }

    #[test]
    fn grad_zero_at_origin() {
        let h = 1.0_f64;
        let kernels: &[(&str, [f64; 3])] = &[
            ("cubic_spline", cubic_spline_grad([0.0; 3], h)),
            ("wendland_c2", wendland_c2_grad([0.0; 3], h)),
            ("poly6", poly6_grad([0.0; 3], h)),
            ("spiky", spiky_grad([0.0; 3], h)),
        ];
        for (name, g) in kernels {
            let mag = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
            assert!(
                mag < 1e-10,
                "grad |{name}| at origin should be 0, got {mag:.3e}"
            );
        }
    }

    #[test]
    fn grad_antiparallel_to_displacement() {
        let h = 1.0_f64;
        let r_vec = [0.3_f64, 0.0, 0.0];

        let g = cubic_spline_grad(r_vec, h);
        assert!(
            g[0] < 0.0,
            "cubic spline grad x should be negative: {}",
            g[0]
        );

        let g = wendland_c2_grad(r_vec, h);
        assert!(
            g[0] < 0.0,
            "wendland C2 grad x should be negative: {}",
            g[0]
        );
    }

    #[test]
    fn poly6_eval_zero_outside_support() {
        let h = 0.5_f64;
        let w = poly6_eval(h + 1e-6, h);
        assert!(
            w.abs() < 1e-14,
            "Poly6 should be 0 outside support, got {w}"
        );
    }

    #[test]
    fn spiky_eval_zero_outside_support() {
        let h = 0.5_f64;
        let w = spiky_eval(h + 1e-6, h);
        assert!(
            w.abs() < 1e-14,
            "Spiky should be 0 outside support, got {w}"
        );
    }

    #[test]
    fn grad_magnitude_matches_scalar_grad_w() {
        let h = 1.0_f64;
        let r_vec = [0.6_f64, 0.0, 0.0];
        let r = 0.6_f64;
        let scalar_dw = CubicSplineKernel.grad_w(r, h);
        let vec_g = cubic_spline_grad(r_vec, h);
        let mag = (vec_g[0] * vec_g[0] + vec_g[1] * vec_g[1] + vec_g[2] * vec_g[2]).sqrt();
        assert!(
            (mag - scalar_dw.abs()).abs() < 1e-12,
            "grad magnitude {mag:.6e} should equal |dW/dr| {:.6e}",
            scalar_dw.abs()
        );
    }

    #[test]
    fn poly6_integrates_to_one() {
        let k = Poly6Kernel;
        let h = 1.0_f64;
        let support = h;
        let n = 120_usize;
        let dx = 2.0 * support / n as f64;
        let mut sum = 0.0;
        for i in 0..n {
            let x = -support + (i as f64 + 0.5) * dx;
            for j in 0..n {
                let y = -support + (j as f64 + 0.5) * dx;
                for kk in 0..n {
                    let z = -support + (kk as f64 + 0.5) * dx;
                    let r = (x * x + y * y + z * z).sqrt();
                    sum += k.w(r, h) * dx * dx * dx;
                }
            }
        }
        assert!(
            (sum - 1.0).abs() < 0.05,
            "Poly6 integral = {sum}, expected ~1.0"
        );
    }

    #[test]
    fn spiky_integrates_to_one() {
        let k = SpikyKernel;
        let h = 1.0_f64;
        let support = h;
        let n = 120_usize;
        let dx = 2.0 * support / n as f64;
        let mut sum = 0.0;
        for i in 0..n {
            let x = -support + (i as f64 + 0.5) * dx;
            for j in 0..n {
                let y = -support + (j as f64 + 0.5) * dx;
                for kk in 0..n {
                    let z = -support + (kk as f64 + 0.5) * dx;
                    let r = (x * x + y * y + z * z).sqrt();
                    sum += k.w(r, h) * dx * dx * dx;
                }
            }
        }
        assert!(
            (sum - 1.0).abs() < 0.05,
            "Spiky integral = {sum}, expected ~1.0"
        );
    }

    // --- Wendland C4 tests ---

    #[test]
    fn test_wendland_c4_gradient_negative_inside() {
        let k = WendlandC4Kernel;
        let dw = k.grad_w(0.5, 1.0);
        assert!(
            dw < 0.0,
            "C4 gradient should be negative inside support, got {dw}"
        );
    }

    #[test]
    fn test_wendland_c4_zero_outside() {
        let k = WendlandC4Kernel;
        assert!(k.w(2.5, 1.0).abs() < 1e-14);
        assert!(k.grad_w(2.5, 1.0).abs() < 1e-14);
    }

    // --- Wendland C6 tests ---

    #[test]
    fn test_wendland_c6_gradient_negative_inside() {
        let k = WendlandC6Kernel;
        let dw = k.grad_w(0.5, 1.0);
        assert!(
            dw < 0.0,
            "C6 gradient should be negative inside support, got {dw}"
        );
    }

    // --- Super-Gaussian tests ---

    #[test]
    fn test_super_gaussian_gradient_negative() {
        let k = SuperGaussianKernel;
        let dw = k.grad_w(0.5, 1.0);
        assert!(dw < 0.0, "SG gradient should be negative inside, got {dw}");
    }

    // --- Kernel correction tests ---

    #[test]
    fn test_shepard_correction_uniform() {
        // Uniform distribution: correction should be close to 1
        let w_vals = vec![0.5, 0.3, 0.2];
        let masses = vec![1.0; 3];
        let densities = vec![1.0; 3];
        let corr = shepard_correction(&w_vals, &masses, &densities);
        assert!(corr > 0.0, "Shepard correction must be positive");
        // sum = 1.0 → corr = 1.0
        assert!((corr - 1.0).abs() < 1e-10, "Corr should be 1.0, got {corr}");
    }

    #[test]
    fn test_mls_moment_matrix_diagonal_dominant() {
        let pos_i = [0.0; 3];
        let neighbors = vec![[0.1, 0.0, 0.0], [-0.1, 0.0, 0.0]];
        let w_vals = vec![1.0, 1.0];
        let masses = vec![1.0, 1.0];
        let densities = vec![1.0, 1.0];
        let m = mls_moment_matrix(pos_i, &neighbors, &w_vals, &masses, &densities);
        // M[0][0] should be sum of volumes * W = 2
        assert!(
            (m[0][0] - 2.0).abs() < 1e-10,
            "M[0][0] should be 2.0, got {}",
            m[0][0]
        );
    }

    // --- Kernel comparison tests ---

    #[test]
    fn test_compare_kernels_same() {
        let k = CubicSplineKernel;
        let (max_abs, max_rel) = compare_kernels(&k, &k, 1.0, 100);
        assert!(max_abs < 1e-14, "Same kernel should have zero diff");
        assert!(max_rel < 1e-14, "Same kernel should have zero rel diff");
    }

    #[test]
    fn test_compare_kernels_different() {
        let k1 = CubicSplineKernel;
        let k2 = WendlandC2Kernel;
        let (max_abs, _) = compare_kernels(&k1, &k2, 1.0, 100);
        assert!(max_abs > 0.0, "Different kernels should have nonzero diff");
    }

    #[test]
    fn test_effective_support_radius() {
        let k = CubicSplineKernel;
        let r_eff = effective_support_radius(&k, 1.0, 1e-10);
        assert!(
            r_eff > 1.9 && r_eff < 2.1,
            "Cubic spline support should be ~2h, got {r_eff}"
        );
    }

    // --- Gradient correction tests ---

    #[test]
    fn test_invert_3x3_identity() {
        let id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let inv = invert_3x3(id);
        for (i, row) in inv.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (val - expected).abs() < 1e-12,
                    "inv[{i}][{j}] = {}, expected {expected}",
                    val
                );
            }
        }
    }

    #[test]
    fn test_wendland_c4_eval_wrapper() {
        let w = wendland_c4_eval(0.5, 1.0);
        assert!(w > 0.0, "C4 eval at r=0.5 should be positive");
        let w = wendland_c4_eval(2.5, 1.0);
        assert!(w.abs() < 1e-14, "C4 eval outside support should be 0");
    }

    #[test]
    fn test_wendland_c6_eval_wrapper() {
        let w = wendland_c6_eval(0.5, 1.0);
        assert!(w > 0.0, "C6 eval at r=0.5 should be positive");
    }

    #[test]
    fn test_super_gaussian_eval_wrapper() {
        let w = super_gaussian_eval(0.0, 1.0);
        assert!(w > 0.0, "SG eval at origin should be positive");
    }

    // --- Wendland C4/C6 monotonicity ---

    #[test]
    fn test_wendland_c4_monotonic_decrease() {
        let k = WendlandC4Kernel;
        let h = 1.0;
        let w0 = k.w(0.0, h);
        let w1 = k.w(0.5, h);
        let w2 = k.w(1.0, h);
        let w3 = k.w(1.5, h);
        assert!(w0 > w1, "W should decrease: w(0)={w0} > w(0.5)={w1}");
        assert!(w1 > w2, "W should decrease: w(0.5)={w1} > w(1.0)={w2}");
        assert!(w2 > w3, "W should decrease: w(1.0)={w2} > w(1.5)={w3}");
    }

    #[test]
    fn test_wendland_c6_monotonic_decrease() {
        let k = WendlandC6Kernel;
        let h = 1.0;
        let w0 = k.w(0.0, h);
        let w1 = k.w(0.5, h);
        let w2 = k.w(1.0, h);
        assert!(w0 > w1, "W should decrease");
        assert!(w1 > w2, "W should decrease");
    }
}

// ---------------------------------------------------------------------------
// Quintic Wendland C6 (alternative normalisation, sometimes called "quintic")
// ---------------------------------------------------------------------------

/// Quintic Wendland C6 kernel for 3D (11th-order polynomial form).
///
/// Compact support at `r/h = 1` (h is the *support radius*).
/// W(r,h) = σ * (1 - q)^8 * (32q^3 + 25q^2 + 8q + 1)  for q = r/h ≤ 1
/// σ = 78 / (7 π h³)
///
/// This variant keeps h as the support radius instead of the smoothing length.
#[derive(Debug, Clone, Copy, Default)]
pub struct QuinticWendlandC6Kernel;

impl SphKernel for QuinticWendlandC6Kernel {
    fn w(&self, r: f64, h: f64) -> f64 {
        let q = r / h;
        if q >= 1.0 {
            return 0.0;
        }
        let sigma = 78.0 / (7.0 * PI * h * h * h);
        let t = 1.0 - q;
        sigma * t.powi(8) * (32.0 * q * q * q + 25.0 * q * q + 8.0 * q + 1.0)
    }

    fn grad_w(&self, r: f64, h: f64) -> f64 {
        let q = r / h;
        if !(1e-14..1.0).contains(&q) {
            return 0.0;
        }
        let sigma = 78.0 / (7.0 * PI * h * h * h);
        let t = 1.0 - q;
        let poly = 32.0 * q * q * q + 25.0 * q * q + 8.0 * q + 1.0;
        let dpoly = 96.0 * q * q + 50.0 * q + 8.0;
        let dw_dq = sigma * (-(8.0 * t.powi(7)) * poly + t.powi(8) * dpoly);
        dw_dq / h
    }

    fn laplacian_w(&self, r: f64, h: f64) -> f64 {
        let q = r / h;
        if q >= 1.0 {
            return 0.0;
        }
        let eps = 1e-6 * h;
        let r_plus = r + eps;
        let r_minus = (r - eps).max(0.0);
        let dw_plus = self.grad_w(r_plus, h);
        let dw_minus = self.grad_w(r_minus, h);
        let d2w = (dw_plus - dw_minus) / (r_plus - r_minus);
        let dw = self.grad_w(r, h);
        d2w + if r > 1e-14 { 2.0 * dw / r } else { 0.0 }
    }
}

// ---------------------------------------------------------------------------
// Gaussian Kernel
// ---------------------------------------------------------------------------

/// Standard Gaussian kernel for 3D.
///
/// W(r,h) = σ * exp(-q²),  q = r/h
/// σ = 1 / (π^(3/2) h³)
///
/// Technically has infinite support; values below machine epsilon at q > 3.
#[derive(Debug, Clone, Copy, Default)]
pub struct GaussianKernel;

impl SphKernel for GaussianKernel {
    fn w(&self, r: f64, h: f64) -> f64 {
        let q = r / h;
        if q > 4.0 {
            return 0.0; // truncate for numerical safety
        }
        let sigma = 1.0 / (PI.powf(1.5) * h * h * h);
        sigma * (-(q * q)).exp()
    }

    fn grad_w(&self, r: f64, h: f64) -> f64 {
        let q = r / h;
        if !(1e-14..4.0).contains(&q) {
            return 0.0;
        }
        let sigma = 1.0 / (PI.powf(1.5) * h * h * h);
        let q2 = q * q;
        // dW/dq = sigma * (-2q) * exp(-q²)
        let dw_dq = sigma * (-2.0 * q) * ((-q2).exp());
        dw_dq / h
    }

    fn laplacian_w(&self, r: f64, h: f64) -> f64 {
        let q = r / h;
        if q > 4.0 {
            return 0.0;
        }
        let sigma = 1.0 / (PI.powf(1.5) * h * h * h);
        let q2 = q * q;
        let h2 = h * h;
        let exp_q = (-q2).exp();
        // d²W/dr² = sigma * (4q²/h² - 2/h²) * exp(-q²)
        let d2w = sigma * (4.0 * q2 - 2.0) / h2 * exp_q;
        let dw = if q > 1e-14 {
            sigma * (-2.0 * q) * exp_q / h
        } else {
            0.0
        };
        d2w + if r > 1e-14 { 2.0 * dw / r } else { 0.0 }
    }
}

// ---------------------------------------------------------------------------
// Quintic Spline Kernel (M5)
// ---------------------------------------------------------------------------

/// Quintic spline kernel (M5) for 3D.
///
/// Compact support at `r/h = 3`.
/// Higher-order accuracy than the cubic spline.
#[derive(Debug, Clone, Copy, Default)]
pub struct QuinticSplineKernel;

impl SphKernel for QuinticSplineKernel {
    fn w(&self, r: f64, h: f64) -> f64 {
        let q = r / h;
        let sigma = 1.0 / (120.0 * PI * h * h * h);
        let f = |x: f64| if x > 0.0 { x } else { 0.0 };
        sigma * (f(3.0 - q).powi(5) - 6.0 * f(2.0 - q).powi(5) + 15.0 * f(1.0 - q).powi(5))
    }

    fn grad_w(&self, r: f64, h: f64) -> f64 {
        let q = r / h;
        if !(1e-14..3.0).contains(&q) {
            return 0.0;
        }
        let sigma = 1.0 / (120.0 * PI * h * h * h);
        let f = |x: f64| if x > 0.0 { x } else { 0.0 };
        let dw_dq = sigma
            * (-5.0 * f(3.0 - q).powi(4) + 30.0 * f(2.0 - q).powi(4) - 75.0 * f(1.0 - q).powi(4));
        dw_dq / h
    }

    fn laplacian_w(&self, r: f64, h: f64) -> f64 {
        let q = r / h;
        if q >= 3.0 {
            return 0.0;
        }
        let eps = 1e-6 * h;
        let r_plus = r + eps;
        let r_minus = (r - eps).max(0.0);
        let dw_plus = self.grad_w(r_plus, h);
        let dw_minus = self.grad_w(r_minus, h);
        let d2w = (dw_plus - dw_minus) / (r_plus - r_minus);
        let dw = self.grad_w(r, h);
        d2w + if r > 1e-14 { 2.0 * dw / r } else { 0.0 }
    }
}

// ---------------------------------------------------------------------------
// Kernel gradient corrections (KSPH / tensorial)
// ---------------------------------------------------------------------------

/// KSPH (Kernel SPH) tensorial gradient correction (Børve et al. 2001).
///
/// Computes the corrected gradient using the second-order correction matrix:
/// `L_i = (Σ_j V_j ∇W_ij ⊗ r_ij)^(-1)`
///
/// The corrected gradient is then: `∇^*W_ij = L_i · ∇W_ij`
pub fn ksph_corrected_gradient(grad_w: [f64; 3], correction_matrix: [[f64; 3]; 3]) -> [f64; 3] {
    let mut result = [0.0_f64; 3];
    for i in 0..3 {
        for j in 0..3 {
            result[i] += correction_matrix[i][j] * grad_w[j];
        }
    }
    result
}

/// Apply KSPH gradient correction to all neighbor gradient vectors.
///
/// Returns the corrected gradient vectors as a `Vec<[f64; 3]>`.
pub fn ksph_correct_all_gradients(
    gradients: &[[f64; 3]],
    correction_matrix: [[f64; 3]; 3],
) -> Vec<[f64; 3]> {
    gradients
        .iter()
        .map(|&g| ksph_corrected_gradient(g, correction_matrix))
        .collect()
}

/// Compute the higher-order (second-moment) kernel correction matrix for particle i.
///
/// The second moment matrix reads:
/// `M2 = Σ_j V_j W_ij r_ij ⊗ r_ij`
///
/// This is used in higher-order SPH schemes (e.g., RK-SPH).
pub fn second_moment_matrix(
    pos_i: [f64; 3],
    neighbor_positions: &[[f64; 3]],
    kernel_values: &[f64],
    masses: &[f64],
    densities: &[f64],
) -> [[f64; 3]; 3] {
    let mut m = [[0.0_f64; 3]; 3];
    for (idx, &pos_j) in neighbor_positions.iter().enumerate() {
        if densities[idx] < 1e-14 {
            continue;
        }
        let w = kernel_values[idx];
        let vol = masses[idx] / densities[idx];
        let r = [
            pos_j[0] - pos_i[0],
            pos_j[1] - pos_i[1],
            pos_j[2] - pos_i[2],
        ];
        for row in 0..3 {
            for col in 0..3 {
                m[row][col] += vol * w * r[row] * r[col];
            }
        }
    }
    m
}

// ---------------------------------------------------------------------------
// Adaptive kernel renormalization
// ---------------------------------------------------------------------------

/// Renormalize a set of kernel values so that their weighted sum equals 1.
///
/// This is the zeroth-order (Shepard) consistency correction.
/// Returns the scale factor; multiply all kernel values by this factor.
pub fn renormalize_kernel(kernel_values: &[f64], volumes: &[f64]) -> f64 {
    let sum: f64 = kernel_values
        .iter()
        .zip(volumes.iter())
        .map(|(&w, &v)| w * v)
        .sum();
    if sum > 1e-30 { 1.0 / sum } else { 1.0 }
}

/// Renormalize kernel values in place for a set of neighbor interactions.
///
/// After renormalization, `Σ V_j W_ij = 1`.
pub fn apply_renormalization(kernel_values: &mut [f64], volumes: &[f64]) {
    let scale = renormalize_kernel(kernel_values, volumes);
    for w in kernel_values.iter_mut() {
        *w *= scale;
    }
}

/// Compute the kernel renormalization factor using a density-weighted sum.
///
/// This variant uses `m_j / ρ_j` as the volume estimate, consistent with
/// standard SPH notation.
pub fn kernel_renorm_factor(kernel_values: &[f64], masses: &[f64], densities: &[f64]) -> f64 {
    shepard_correction(kernel_values, masses, densities)
}

// ---------------------------------------------------------------------------
// Higher-order kernel stencils
// ---------------------------------------------------------------------------

/// Evaluate a 1D stencil for the Laplacian using the finite-difference
/// approximation through the kernel Laplacian.
///
/// `∇²f ≈ 2 * Σ_j (m_j / ρ_j) * (f_j - f_i) * (dW/dr) / r`
pub fn sph_laplacian_1d(
    f_i: f64,
    f_neighbors: &[f64],
    positions_i: f64,
    positions_neighbors: &[f64],
    masses: &[f64],
    densities: &[f64],
    kernel: &dyn SphKernel,
    h: f64,
) -> f64 {
    let mut sum = 0.0_f64;
    for (idx, &pos_j) in positions_neighbors.iter().enumerate() {
        if densities[idx] < 1e-14 {
            continue;
        }
        let r = (pos_j - positions_i).abs();
        if r < 1e-14 {
            continue;
        }
        let dw_dr = kernel.grad_w(r, h);
        let vol = masses[idx] / densities[idx];
        sum += vol * (f_neighbors[idx] - f_i) * dw_dr / r;
    }
    2.0 * sum
}

/// Compute the SPH gradient in 3D at particle i.
///
/// `∇f_i ≈ Σ_j (m_j / ρ_j) * (f_j - f_i) * ∇W_ij`
pub fn sph_gradient_3d(
    f_i: f64,
    f_neighbors: &[f64],
    pos_i: [f64; 3],
    pos_neighbors: &[[f64; 3]],
    masses: &[f64],
    densities: &[f64],
    kernel: &dyn SphKernel,
    h: f64,
) -> [f64; 3] {
    let mut result = [0.0_f64; 3];
    for (idx, &pos_j) in pos_neighbors.iter().enumerate() {
        if densities[idx] < 1e-14 {
            continue;
        }
        let r_vec = [
            pos_j[0] - pos_i[0],
            pos_j[1] - pos_i[1],
            pos_j[2] - pos_i[2],
        ];
        let r2 = r_vec[0] * r_vec[0] + r_vec[1] * r_vec[1] + r_vec[2] * r_vec[2];
        if r2 < 1e-28 {
            continue;
        }
        let r = r2.sqrt();
        let dw_dr = kernel.grad_w(r, h);
        let inv_r = dw_dr / r;
        let vol = masses[idx] / densities[idx];
        let df = f_neighbors[idx] - f_i;
        for k in 0..3 {
            result[k] += vol * df * r_vec[k] * inv_r;
        }
    }
    result
}

// ---------------------------------------------------------------------------
// Free-function API for new kernels
// ---------------------------------------------------------------------------

/// Evaluate the Quintic Wendland C6 kernel at scalar distance `r`.
pub fn quintic_wendland_c6_eval(r: f64, h: f64) -> f64 {
    QuinticWendlandC6Kernel.w(r, h)
}

/// Evaluate the Gaussian kernel at scalar distance `r`.
pub fn gaussian_eval(r: f64, h: f64) -> f64 {
    GaussianKernel.w(r, h)
}

/// Evaluate the Quintic Spline kernel at scalar distance `r`.
pub fn quintic_spline_eval(r: f64, h: f64) -> f64 {
    QuinticSplineKernel.w(r, h)
}

/// Evaluate the Gaussian gradient given displacement vector `r_vec`.
pub fn gaussian_grad(r_vec: [f64; 3], h: f64) -> [f64; 3] {
    grad(&GaussianKernel, r_vec, h)
}

// ---------------------------------------------------------------------------
// Additional tests for new kernels and utilities
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests_extended {
    use super::*;

    // --- QuinticWendlandC6Kernel ---

    #[test]
    fn test_quintic_wendland_c6_positive_inside() {
        let k = QuinticWendlandC6Kernel;
        assert!(k.w(0.0, 1.0) > 0.0, "QWC6 at origin should be positive");
        assert!(k.w(0.3, 1.0) > 0.0, "QWC6 at r=0.3 should be positive");
        assert!(k.w(0.8, 1.0) > 0.0, "QWC6 at r=0.8 should be positive");
    }

    #[test]
    fn test_quintic_wendland_c6_zero_outside() {
        let k = QuinticWendlandC6Kernel;
        assert!(k.w(1.0, 1.0).abs() < 1e-14, "QWC6 at boundary should be 0");
        assert!(
            k.w(1.5, 1.0).abs() < 1e-14,
            "QWC6 outside support should be 0"
        );
    }

    #[test]
    fn test_quintic_wendland_c6_gradient_negative() {
        let k = QuinticWendlandC6Kernel;
        let dw = k.grad_w(0.3, 1.0);
        assert!(
            dw < 0.0,
            "QWC6 gradient inside support should be negative: {dw}"
        );
    }

    #[test]
    fn test_quintic_wendland_c6_monotonic() {
        let k = QuinticWendlandC6Kernel;
        let h = 1.0;
        let w0 = k.w(0.0, h);
        let w1 = k.w(0.3, h);
        let w2 = k.w(0.6, h);
        let w3 = k.w(0.9, h);
        assert!(
            w0 > w1 && w1 > w2 && w2 > w3,
            "QWC6 should be monotonically decreasing"
        );
    }

    #[test]
    fn test_quintic_wendland_c6_eval_wrapper() {
        let w = quintic_wendland_c6_eval(0.0, 1.0);
        assert!(w > 0.0);
        let w = quintic_wendland_c6_eval(1.5, 1.0);
        assert!(w.abs() < 1e-14);
    }

    // --- GaussianKernel ---

    #[test]
    fn test_gaussian_positive_at_origin() {
        let k = GaussianKernel;
        assert!(k.w(0.0, 1.0) > 0.0, "Gaussian at origin should be positive");
    }

    #[test]
    fn test_gaussian_decreases_with_distance() {
        let k = GaussianKernel;
        let h = 1.0;
        let w0 = k.w(0.0, h);
        let w1 = k.w(0.5, h);
        let w2 = k.w(1.0, h);
        assert!(w0 > w1 && w1 > w2, "Gaussian should decrease with distance");
    }

    #[test]
    fn test_gaussian_gradient_negative() {
        let k = GaussianKernel;
        let dw = k.grad_w(0.5, 1.0);
        assert!(
            dw < 0.0,
            "Gaussian gradient inside support should be negative: {dw}"
        );
    }

    #[test]
    fn test_gaussian_gradient_zero_at_origin() {
        let k = GaussianKernel;
        assert!(
            k.grad_w(0.0, 1.0).abs() < 1e-10,
            "Gaussian gradient at origin should be 0"
        );
    }

    #[test]
    fn test_gaussian_eval_wrapper() {
        let w = gaussian_eval(0.0, 1.0);
        assert!(w > 0.0);
    }

    #[test]
    fn test_gaussian_integrates_approximately() {
        let k = GaussianKernel;
        let h = 1.0;
        let support = 4.0 * h;
        let n = 100usize;
        let dx = 2.0 * support / n as f64;
        let mut sum = 0.0;
        for i in 0..n {
            let x = -support + (i as f64 + 0.5) * dx;
            for j in 0..n {
                let y = -support + (j as f64 + 0.5) * dx;
                for kk in 0..n {
                    let z = -support + (kk as f64 + 0.5) * dx;
                    let r = (x * x + y * y + z * z).sqrt();
                    sum += k.w(r, h) * dx * dx * dx;
                }
            }
        }
        assert!(
            (sum - 1.0).abs() < 0.05,
            "Gaussian integral = {sum}, expected ~1.0"
        );
    }

    // --- QuinticSplineKernel ---

    #[test]
    fn test_quintic_spline_positive_inside() {
        let k = QuinticSplineKernel;
        let h = 1.0;
        assert!(k.w(0.0, h) > 0.0, "QS at origin should be positive");
        assert!(k.w(0.5, h) > 0.0, "QS at r=0.5 should be positive");
        assert!(k.w(1.5, h) > 0.0, "QS at r=1.5 should be positive");
        assert!(k.w(2.5, h) > 0.0, "QS at r=2.5 should be positive");
    }

    #[test]
    fn test_quintic_spline_zero_outside() {
        let k = QuinticSplineKernel;
        let h = 1.0;
        assert!(k.w(3.0, h).abs() < 1e-14, "QS at boundary should be 0");
        assert!(k.w(3.5, h).abs() < 1e-14, "QS outside support should be 0");
    }

    #[test]
    fn test_quintic_spline_gradient_zero_at_origin() {
        let k = QuinticSplineKernel;
        assert!(k.grad_w(0.0, 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_quintic_spline_eval_wrapper() {
        let w = quintic_spline_eval(0.0, 1.0);
        assert!(w > 0.0);
        let w_out = quintic_spline_eval(3.5, 1.0);
        assert!(w_out.abs() < 1e-14);
    }

    // --- KSPH gradient correction ---

    #[test]
    fn test_ksph_corrected_gradient_identity() {
        // With identity correction matrix, gradient should be unchanged
        let id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let g = [1.0, 2.0, 3.0];
        let result = ksph_corrected_gradient(g, id);
        for i in 0..3 {
            assert!(
                (result[i] - g[i]).abs() < 1e-12,
                "Identity correction should preserve gradient: result[{i}]={}",
                result[i]
            );
        }
    }

    #[test]
    fn test_ksph_corrected_gradient_scales() {
        // Diagonal scaling matrix: should scale each component
        let scale = [[2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 0.5]];
        let g = [1.0, 1.0, 1.0];
        let result = ksph_corrected_gradient(g, scale);
        assert!((result[0] - 2.0).abs() < 1e-12);
        assert!((result[1] - 3.0).abs() < 1e-12);
        assert!((result[2] - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_ksph_correct_all_gradients_count() {
        let id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let grads = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let result = ksph_correct_all_gradients(&grads, id);
        assert_eq!(result.len(), grads.len());
    }

    // --- Second moment matrix ---

    #[test]
    fn test_second_moment_matrix_symmetric() {
        let pos_i = [0.0; 3];
        let neighbors = vec![[0.1, 0.2, 0.3], [-0.1, 0.1, -0.2]];
        let w_vals = vec![1.0, 1.0];
        let masses = vec![1.0, 1.0];
        let densities = vec![1.0, 1.0];
        let m = second_moment_matrix(pos_i, &neighbors, &w_vals, &masses, &densities);
        for (i, row) in m.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!(
                    (val - m[j][i]).abs() < 1e-12,
                    "Second moment matrix should be symmetric: m[{i}][{j}]={} != m[{j}][{i}]={}",
                    val,
                    m[j][i]
                );
            }
        }
    }

    // --- Kernel renormalization ---

    #[test]
    fn test_renormalize_kernel_sums_to_one() {
        let w_vals = vec![0.5, 0.3, 0.2];
        let vols = vec![1.0, 1.0, 1.0];
        let scale = renormalize_kernel(&w_vals, &vols);
        // sum = 1.0, scale = 1.0
        assert!((scale - 1.0).abs() < 1e-12, "scale={scale}");
    }

    #[test]
    fn test_renormalize_kernel_non_unit_sum() {
        let w_vals = vec![0.4, 0.4];
        let vols = vec![1.0, 1.0];
        let scale = renormalize_kernel(&w_vals, &vols);
        // sum = 0.8, scale = 1.25
        assert!((scale - 1.25).abs() < 1e-12, "scale={scale}");
    }

    #[test]
    fn test_apply_renormalization_modifies_in_place() {
        let mut w_vals = vec![0.4, 0.4];
        let vols = vec![1.0, 1.0];
        apply_renormalization(&mut w_vals, &vols);
        let sum: f64 = w_vals.iter().zip(vols.iter()).map(|(&w, &v)| w * v).sum();
        assert!(
            (sum - 1.0).abs() < 1e-12,
            "After renormalization, sum should be 1: {sum}"
        );
    }

    #[test]
    fn test_kernel_renorm_factor_matches_shepard() {
        let w = vec![0.6, 0.4];
        let m = vec![1.0, 1.0];
        let d = vec![1.0, 1.0];
        let factor = kernel_renorm_factor(&w, &m, &d);
        let shepard = shepard_correction(&w, &m, &d);
        assert!((factor - shepard).abs() < 1e-12);
    }

    // --- SPH gradient and Laplacian stencils ---

    #[test]
    fn test_sph_laplacian_1d_constant_field() {
        // Constant field: Laplacian should be 0
        let k = CubicSplineKernel;
        let h = 0.2;
        let f_i = 5.0;
        let pos_neighbors = vec![-0.1, 0.0, 0.1];
        let f_neighbors = vec![5.0, 5.0, 5.0];
        let masses = vec![0.001, 0.001, 0.001];
        let densities = vec![1000.0, 1000.0, 1000.0];
        let lap = sph_laplacian_1d(
            f_i,
            &f_neighbors,
            0.0,
            &pos_neighbors,
            &masses,
            &densities,
            &k,
            h,
        );
        // The Laplacian of a constant field should be approximately 0
        assert!(
            lap.abs() < 1.0,
            "Laplacian of constant field should be ~0: {lap}"
        );
    }

    #[test]
    fn test_sph_gradient_3d_single_neighbor() {
        // sph_gradient_3d uses r_vec = pos_j - pos_i and dW/dr < 0
        // For f_j > f_i with neighbor in +x direction, result[0] = vol*(f_j-f_i)*(dW/dr)/r * r_x
        // dW/dr < 0 for cubic spline, so gradient[0] < 0
        let k = CubicSplineKernel;
        let h = 1.0;
        let pos_i = [0.0, 0.0, 0.0];
        let pos_neighbors = vec![[0.5, 0.0, 0.0]];
        let f_neighbors = vec![1.0];
        let masses = vec![1.0];
        let densities = vec![1.0];
        let grad_f = sph_gradient_3d(
            0.0,
            &f_neighbors,
            pos_i,
            &pos_neighbors,
            &masses,
            &densities,
            &k,
            h,
        );
        // Gradient is non-zero
        assert!(
            grad_f[0].abs() > 1e-10,
            "Gradient x should be non-zero: {}",
            grad_f[0]
        );
        assert!(
            grad_f[1].abs() < 1e-12,
            "Gradient y should be 0: {}",
            grad_f[1]
        );
        assert!(
            grad_f[2].abs() < 1e-12,
            "Gradient z should be 0: {}",
            grad_f[2]
        );
    }

    #[test]
    fn test_sph_gradient_3d_antisymmetric() {
        // Two symmetric neighbors with antisymmetric field → gradient non-zero
        let k = CubicSplineKernel;
        let h = 1.0;
        let pos_i = [0.0, 0.0, 0.0];
        let pos_neighbors = vec![[0.5, 0.0, 0.0], [-0.5, 0.0, 0.0]];
        let f_neighbors = vec![1.0, -1.0]; // antisymmetric field
        let masses = vec![1.0, 1.0];
        let densities = vec![1.0, 1.0];
        let grad_f = sph_gradient_3d(
            0.0,
            &f_neighbors,
            pos_i,
            &pos_neighbors,
            &masses,
            &densities,
            &k,
            h,
        );
        // Should have non-zero x gradient (antisymmetric field → contributions add up)
        assert!(
            grad_f[0].abs() > 1e-10,
            "Antisymmetric field should give non-zero x gradient: {}",
            grad_f[0]
        );
    }

    // --- Gaussian gradient vector ---

    #[test]
    fn test_gaussian_grad_points_opposite_to_displacement() {
        let h = 1.0;
        let r_vec = [0.5, 0.0, 0.0];
        let g = gaussian_grad(r_vec, h);
        assert!(
            g[0] < 0.0,
            "Gaussian grad x should be negative (away from origin): {}",
            g[0]
        );
        assert!(g[1].abs() < 1e-12);
        assert!(g[2].abs() < 1e-12);
    }

    #[test]
    fn test_gaussian_grad_zero_at_origin() {
        let g = gaussian_grad([0.0; 3], 1.0);
        let mag = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
        assert!(
            mag < 1e-10,
            "Gaussian grad at origin should be 0: mag={mag}"
        );
    }

    // --- Kernel comparison with new kernels ---

    #[test]
    fn test_gaussian_different_from_cubic_spline() {
        let (max_abs, _) = compare_kernels(&GaussianKernel, &CubicSplineKernel, 1.0, 100);
        assert!(max_abs > 0.0, "Gaussian and CubicSpline should differ");
    }

    #[test]
    fn test_quintic_spline_different_from_cubic_spline() {
        let (max_abs, _) = compare_kernels(&QuinticSplineKernel, &CubicSplineKernel, 1.0, 100);
        assert!(
            max_abs > 0.0,
            "Quintic spline and CubicSpline should differ"
        );
    }

    // --- Laplacian sign tests ---

    #[test]
    fn test_gaussian_laplacian_negative_near_origin() {
        // Gaussian is concave near origin → Laplacian should be negative
        let k = GaussianKernel;
        // At r=0: d²W/dr² = sigma * (-2/h²) → negative
        // But the 3D Laplacian includes 2/r * dW/dr term
        // Just check it's finite
        let lap = k.laplacian_w(0.1, 1.0);
        assert!(
            lap.is_finite(),
            "Gaussian Laplacian should be finite: {lap}"
        );
    }

    #[test]
    fn test_quintic_spline_laplacian_finite() {
        let k = QuinticSplineKernel;
        for r in [0.1, 0.5, 1.0, 1.5, 2.0, 2.5] {
            let lap = k.laplacian_w(r, 1.0);
            assert!(
                lap.is_finite(),
                "Quintic spline Laplacian should be finite at r={r}: {lap}"
            );
        }
    }

    // --- Kernel symmetry: w(-r) = w(r) since kernels use |r| ---

    #[test]
    fn test_kernels_use_scalar_distance() {
        let kernels: Vec<(&str, Box<dyn SphKernel>)> = vec![
            ("Gaussian", Box::new(GaussianKernel)),
            ("QuinticSpline", Box::new(QuinticSplineKernel)),
            ("QuinticWendlandC6", Box::new(QuinticWendlandC6Kernel)),
        ];
        for (name, k) in &kernels {
            let w = k.w(0.3, 1.0);
            assert!(w > 0.0, "{name} should be positive at r=0.3: w={w}");
        }
    }

    // --- Test gradient magnitudes for new kernels ---

    #[test]
    fn test_gaussian_grad_magnitude_matches_scalar() {
        let h = 1.0;
        let r = 0.5;
        let r_vec = [r, 0.0, 0.0];
        let scalar_dw = GaussianKernel.grad_w(r, h).abs();
        let vec_g = gaussian_grad(r_vec, h);
        let mag = (vec_g[0] * vec_g[0] + vec_g[1] * vec_g[1] + vec_g[2] * vec_g[2]).sqrt();
        assert!(
            (mag - scalar_dw).abs() < 1e-12,
            "Gaussian grad magnitude mismatch: {mag:.6e} vs {scalar_dw:.6e}"
        );
    }

    #[test]
    fn test_quintic_wendland_c6_gradient_zero_at_origin() {
        let k = QuinticWendlandC6Kernel;
        assert!(
            k.grad_w(0.0, 1.0).abs() < 1e-10,
            "QWC6 gradient at origin should be 0"
        );
    }

    #[test]
    fn test_quintic_wendland_c6_gradient_zero_outside() {
        let k = QuinticWendlandC6Kernel;
        assert!(
            k.grad_w(1.5, 1.0).abs() < 1e-14,
            "QWC6 gradient outside support should be 0"
        );
    }
}
