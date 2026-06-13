// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Chebyshev polynomial smoother for AMG cycles.
//!
//! A Chebyshev smoother damps the high-frequency (large-eigenvalue) part of the
//! error spectrum with a fixed-degree polynomial built from a spectral-interval
//! estimate.  Unlike Gauss-Seidel it has no sequential row dependency, so it
//! parallelises trivially, and its smoothing behaviour is predictable and
//! mesh-independent because the polynomial is constructed analytically.
//!
//! ## Operator
//!
//! The smoother is *Jacobi-preconditioned*: it builds a polynomial in the
//! diagonally-scaled operator `D⁻¹A`.  This is the standard AMG Chebyshev
//! smoother (Adams, Brezina, Hu & Tuminaro 2003) and is exactly the operator
//! whose spectral radius [`power_iteration_spectral_radius`] estimates, so the
//! eigenvalue bound feeds the recurrence directly.
//!
//! ## Spectral interval
//!
//! The polynomial is optimised over the band `[λmax / θ, λmax]` with `θ ≈ 30`,
//! so it targets the top `1/θ` of the spectrum — precisely the components a
//! coarse grid cannot represent.  `λmax` is the estimated spectral radius of
//! `D⁻¹A`, inflated by [`SPECTRAL_SAFETY`] so the whole upper spectrum lies
//! inside the band (a slight over-estimate cannot amplify any mode).
//!
//! ## Recurrence (3-term)
//!
//! With `θ = (λmax + λmin)/2`, `δ = (λmax − λmin)/2`, `σ = θ/δ`, `ρ₀ = 1/σ`:
//!
//! ```text
//! r  = b − A·x
//! d  = D⁻¹r / θ
//! repeat `degree` times:
//!     x      += d
//!     r       = b − A·x
//!     ρ_new   = 1 / (2σ − ρ)
//!     d       = (ρ·ρ_new)·d + (2·ρ_new/δ)·D⁻¹r
//!     ρ       = ρ_new
//! ```
//!
//! After `degree` iterations the error is `p(D⁻¹A)·e₀` with `deg p = degree`,
//! `p(0) = 1`, minimising `maxₓ∈band |p(x)|` — the classical min-max property.

use crate::parallel_solver::CsrMatrix;
use crate::solvers::amg::smoothed_aggregation::power_iteration_spectral_radius;

/// Ratio `θ` between the largest and smallest eigenvalue bounds.
///
/// The polynomial is optimised over `[λmax / CHEBYSHEV_THETA, λmax]`, so the
/// smoother concentrates on the top `1/θ` of the spectrum.
const CHEBYSHEV_THETA: f64 = 30.0;

/// Number of power iterations used to estimate `λmax = ρ(D⁻¹A)`.
const SPECTRAL_RADIUS_ITERS: usize = 20;

/// Safety inflation applied to the estimated spectral radius.
///
/// Using `1.1·ρ` guarantees the true `λmax` is contained in the band; a slight
/// over-estimate keeps every mode inside the optimisation interval (where the
/// polynomial is bounded by `1`), whereas an under-estimate could *amplify* the
/// modes above the assumed `λmax`.
const SPECTRAL_SAFETY: f64 = 1.1;

/// Threshold below which a value is treated as numerically zero.
const TINY: f64 = 1e-300;

/// Apply a Jacobi-preconditioned Chebyshev smoother to `A·x = b`.
///
/// `a` is a dense, row-major, square operator (the same role the Gauss-Seidel
/// smoothers play for sparse operators).  `x` is updated in place from its
/// current value (use a zero vector for a cold start).  `degree` is the
/// polynomial degree, i.e. the number of smoothing sweeps.
///
/// The spectral radius of `D⁻¹A` is estimated internally via
/// [`power_iteration_spectral_radius`]; `λmin` is set to `λmax / θ`.
pub fn chebyshev_smoother(a: &[Vec<f64>], b: &[f64], x: &mut [f64], degree: usize) {
    let n = a.len();
    if n == 0 || degree == 0 {
        return;
    }

    // Estimate λmax = ρ(D⁻¹A) (inflated for safety) and form the band.
    let csr = dense_to_csr(a);
    let lambda_max = SPECTRAL_SAFETY * power_iteration_spectral_radius(&csr, SPECTRAL_RADIUS_ITERS);
    let lambda_min = lambda_max / CHEBYSHEV_THETA;

    let theta = 0.5 * (lambda_max + lambda_min);
    let delta = 0.5 * (lambda_max - lambda_min);
    if delta.abs() < TINY || theta.abs() < TINY {
        return;
    }
    let sigma = theta / delta;
    let mut rho = 1.0 / sigma;

    // Reciprocal diagonal D⁻¹ (Jacobi preconditioner).
    let d_inv: Vec<f64> = a
        .iter()
        .enumerate()
        .map(|(i, row)| {
            let d = row[i];
            if d.abs() > TINY { 1.0 / d } else { 1.0 }
        })
        .collect();

    // r = b − A·x ; d = D⁻¹r / θ.
    let mut ax = vec![0.0f64; n];
    dense_matvec(a, x, &mut ax);
    let mut dir: Vec<f64> = d_inv
        .iter()
        .zip(b.iter())
        .zip(ax.iter())
        .map(|((di, bi), axi)| di * (bi - axi) / theta)
        .collect();

    for _ in 0..degree {
        // x += d.
        for (xi, di) in x.iter_mut().zip(dir.iter()) {
            *xi += di;
        }

        // r = b − A·x.
        dense_matvec(a, x, &mut ax);

        let rho_new = 1.0 / (2.0 * sigma - rho);
        let coef_dir = rho * rho_new;
        let coef_res = 2.0 * rho_new / delta;

        // d = coef_dir·d + coef_res·D⁻¹r.
        for ((di, di_inv), (bi, axi)) in dir
            .iter_mut()
            .zip(d_inv.iter())
            .zip(b.iter().zip(ax.iter()))
        {
            let z = di_inv * (bi - axi);
            *di = coef_dir * *di + coef_res * z;
        }
        rho = rho_new;
    }
}

/// Dense row-major matrix-vector product `y = A·x`.
fn dense_matvec(a: &[Vec<f64>], x: &[f64], y: &mut [f64]) {
    for (yi, row) in y.iter_mut().zip(a.iter()) {
        *yi = row.iter().zip(x.iter()).map(|(aij, xj)| aij * xj).sum();
    }
}

/// Convert a dense row-major matrix to CSR, dropping exact zeros.
fn dense_to_csr(a: &[Vec<f64>]) -> CsrMatrix {
    let n = a.len();
    let ncols = a.first().map_or(0, Vec::len);
    let mut row_offsets = vec![0usize; n + 1];
    let mut col_indices = Vec::new();
    let mut values = Vec::new();
    for (i, row) in a.iter().enumerate() {
        for (j, &v) in row.iter().enumerate() {
            if v != 0.0 {
                col_indices.push(j);
                values.push(v);
            }
        }
        row_offsets[i + 1] = col_indices.len();
    }
    CsrMatrix {
        nrows: n,
        ncols,
        row_offsets,
        col_indices,
        values,
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Dense 1-D Poisson tridiagonal `[-1, 2, -1]` of size `n`.
    fn poisson_1d_dense(n: usize) -> Vec<Vec<f64>> {
        let mut a = vec![vec![0.0f64; n]; n];
        for (i, row) in a.iter_mut().enumerate() {
            row[i] = 2.0;
            if i > 0 {
                row[i - 1] = -1.0;
            }
            if i + 1 < n {
                row[i + 1] = -1.0;
            }
        }
        a
    }

    fn l2(v: &[f64]) -> f64 {
        v.iter().map(|x| x * x).sum::<f64>().sqrt()
    }

    #[test]
    fn chebyshev_damps_high_frequency() {
        // b = 0 ⇒ exact solution is 0, so x is the error itself.
        let n = 33;
        let a = poisson_1d_dense(n);
        let b = vec![0.0f64; n];
        let x0: Vec<f64> = (0..n).map(|i| if i % 2 == 0 { 1.0 } else { -1.0 }).collect();
        let norm0 = l2(&x0);

        let mut x = x0.clone();
        chebyshev_smoother(&a, &b, &mut x, 5);
        let factor = l2(&x) / norm0;

        assert!(
            factor < 0.6,
            "Chebyshev did not damp the high-frequency mode: factor = {factor:.4}"
        );
    }

    #[test]
    fn chebyshev_reduces_residual() {
        // Generic RHS: the smoother must reduce the residual ‖b − A·x‖.
        let n = 17;
        let a = poisson_1d_dense(n);
        let b: Vec<f64> = (0..n).map(|i| (i as f64 * 0.3).sin()).collect();
        let mut x = vec![0.0f64; n];

        let mut ax = vec![0.0f64; n];
        dense_matvec(&a, &x, &mut ax);
        let r0 = l2(&b.iter().zip(ax.iter()).map(|(bi, ai)| bi - ai).collect::<Vec<_>>());

        chebyshev_smoother(&a, &b, &mut x, 8);

        dense_matvec(&a, &x, &mut ax);
        let r1 = l2(&b.iter().zip(ax.iter()).map(|(bi, ai)| bi - ai).collect::<Vec<_>>());

        assert!(r1 < r0, "Chebyshev did not reduce residual: {r0} -> {r1}");
    }

    #[test]
    fn higher_degree_smooths_more() {
        // More degree ⇒ stronger damping of the checkerboard error.
        let n = 33;
        let a = poisson_1d_dense(n);
        let b = vec![0.0f64; n];
        let x0: Vec<f64> = (0..n).map(|i| if i % 2 == 0 { 1.0 } else { -1.0 }).collect();
        let norm0 = l2(&x0);

        let mut x2 = x0.clone();
        chebyshev_smoother(&a, &b, &mut x2, 2);
        let mut x6 = x0.clone();
        chebyshev_smoother(&a, &b, &mut x6, 6);

        assert!(
            l2(&x6) / norm0 < l2(&x2) / norm0,
            "degree 6 ({:.4}) should out-damp degree 2 ({:.4})",
            l2(&x6) / norm0,
            l2(&x2) / norm0
        );
    }
}
