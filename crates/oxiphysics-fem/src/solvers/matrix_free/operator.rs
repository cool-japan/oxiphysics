// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Matrix-free linear operators and a preconditioned Conjugate Gradient solver.
//!
//! A *matrix-free* operator applies the action `y = A·x` without ever forming
//! the matrix `A`. This is the foundation of the Kronbichler-Kormann (2012)
//! sum-factorization approach: the global stiffness matrix of a high-order
//! finite-element discretization is never assembled; instead each element's
//! contribution is evaluated on the fly via 1D tensor contractions. The result
//! is an operator whose memory footprint is `O(1)` per element (a handful of
//! small 1D matrices) yet whose mat-vec cost is `O(p^{d+1})` per element rather
//! than the `O(p^{2d})` of a dense element matrix.
//!
//! This module provides:
//! - [`MatrixFreeOperator`]: the trait implemented by such operators.
//! - [`MatrixFreePcg`]: a left-preconditioned Conjugate Gradient solver that
//!   drives any [`MatrixFreeOperator`] with any
//!   [`crate::solvers::amg::preconditioner::Preconditioner`].

use crate::solvers::amg::preconditioner::Preconditioner;

/// A matrix-free linear operator: `y = A·x` without forming `A` explicitly.
pub trait MatrixFreeOperator: Send + Sync {
    /// Apply the operator: `y += A·x`. The caller is responsible for zeroing `y`
    /// beforehand if a plain `y = A·x` is wanted.
    fn apply(&self, x: &[f64], y: &mut [f64]);

    /// Dimension of the (square) operator (`n × n`).
    fn dim(&self) -> usize;
}

/// Error type for the matrix-free PCG solver.
#[derive(Debug, thiserror::Error)]
pub enum MatrixFreeError {
    /// CG failed to reach the tolerance within the iteration budget.
    #[error("CG did not converge in {0} iterations (final residual: {1:.3e})")]
    DidNotConverge(usize, f64),
    /// The right-hand side length did not match the operator dimension.
    #[error("dimension mismatch: operator {0}, rhs {1}")]
    DimensionMismatch(usize, usize),
}

/// A matrix-free preconditioned Conjugate Gradient solver.
///
/// Works with any [`MatrixFreeOperator`] and any [`Preconditioner`].
pub struct MatrixFreePcg {
    /// Maximum number of CG iterations.
    pub max_iters: usize,
    /// Relative residual tolerance (‖r‖/‖b‖) for convergence.
    pub tol: f64,
}

/// Euclidean inner product `Σ_i a[i]·b[i]`.
#[inline]
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(ai, bi)| ai * bi).sum()
}

impl MatrixFreePcg {
    /// Create a solver with the given iteration budget and relative tolerance.
    pub fn new(max_iters: usize, tol: f64) -> Self {
        Self { max_iters, tol }
    }

    /// Solve `A·x = b` by left-preconditioned Conjugate Gradient.
    ///
    /// `op` supplies the action of the (symmetric positive-definite) operator
    /// `A`, and `precond` supplies an approximation `M^{-1} ≈ A^{-1}`. On entry
    /// `x` holds the initial guess; on success it holds the solution and the
    /// returned `usize` is the number of iterations taken (0 if the initial
    /// guess already satisfies the tolerance).
    ///
    /// # Errors
    /// Returns [`MatrixFreeError::DimensionMismatch`] if `b` or `x` does not
    /// match `op.dim()`, and [`MatrixFreeError::DidNotConverge`] if the relative
    /// residual stays above `tol` after `max_iters` iterations.
    pub fn solve(
        &self,
        op: &dyn MatrixFreeOperator,
        precond: &dyn Preconditioner,
        b: &[f64],
        x: &mut [f64],
    ) -> Result<usize, MatrixFreeError> {
        let n = op.dim();
        if b.len() != n {
            return Err(MatrixFreeError::DimensionMismatch(n, b.len()));
        }
        if x.len() != n {
            return Err(MatrixFreeError::DimensionMismatch(n, x.len()));
        }

        // Initial residual r = b - A·x. `apply` accumulates, and `ax` starts at
        // zero, so `ax = A·x` after the call.
        let mut ax = vec![0.0f64; n];
        op.apply(x, &mut ax);
        let mut r: Vec<f64> = b.iter().zip(ax.iter()).map(|(bi, axi)| bi - axi).collect();

        // Guard against a zero right-hand side: use a floor so the relative test
        // never divides by zero.
        let b_norm = dot(b, b).sqrt().max(1e-300);

        // Early-out if the initial guess already satisfies the tolerance.
        if dot(&r, &r).sqrt() / b_norm <= self.tol {
            return Ok(0);
        }

        // z = M^{-1} r, search direction p = z.
        let mut z = vec![0.0f64; n];
        precond.apply(&r, &mut z);
        let mut p = z.clone();
        let mut rz_old = dot(&r, &z);

        let mut ap = vec![0.0f64; n];
        for iter in 1..=self.max_iters {
            // ap = A·p (zero first, since `apply` accumulates).
            for ai in ap.iter_mut() {
                *ai = 0.0;
            }
            op.apply(&p, &mut ap);

            let p_ap = dot(&p, &ap);
            if p_ap.abs() < 1e-300 {
                // Breakdown: p is in the null space of A (or numerically so).
                break;
            }
            let alpha = rz_old / p_ap;

            for i in 0..n {
                x[i] += alpha * p[i];
                r[i] -= alpha * ap[i];
            }

            let res = dot(&r, &r).sqrt();
            if res / b_norm <= self.tol {
                return Ok(iter);
            }

            precond.apply(&r, &mut z);
            let rz_new = dot(&r, &z);

            if rz_old.abs() < 1e-300 {
                // Breakdown: cannot form a meaningful conjugacy factor.
                break;
            }
            let beta = rz_new / rz_old;

            for i in 0..n {
                p[i] = z[i] + beta * p[i];
            }
            rz_old = rz_new;
        }

        let res = dot(&r, &r).sqrt();
        Err(MatrixFreeError::DidNotConverge(self.max_iters, res))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A tiny dense operator for testing: `y[i] += Σ_j a[i*n+j]·x[j]`.
    struct DenseOp {
        a: Vec<f64>,
        n: usize,
    }

    impl MatrixFreeOperator for DenseOp {
        fn apply(&self, x: &[f64], y: &mut [f64]) {
            for (i, yi) in y.iter_mut().enumerate().take(self.n) {
                let row = &self.a[i * self.n..i * self.n + self.n];
                let s: f64 = row.iter().zip(x.iter()).map(|(aij, xj)| aij * xj).sum();
                *yi += s;
            }
        }
        fn dim(&self) -> usize {
            self.n
        }
    }

    /// Identity preconditioner: `z = r`.
    struct IdentityPrecond {
        n: usize,
    }

    impl Preconditioner for IdentityPrecond {
        fn apply(&self, r: &[f64], z: &mut [f64]) {
            z.copy_from_slice(r);
        }
        fn n(&self) -> usize {
            self.n
        }
    }

    #[test]
    fn pcg_solves_small_spd_system() {
        // Diagonally dominant SPD matrix:
        //   [[4, 1, 0],
        //    [1, 4, 1],
        //    [0, 1, 3]]
        let op = DenseOp {
            a: vec![4.0, 1.0, 0.0, 1.0, 4.0, 1.0, 0.0, 1.0, 3.0],
            n: 3,
        };
        let precond = IdentityPrecond { n: 3 };
        let b = [1.0, 2.0, 3.0];
        let mut x = [0.0; 3];

        let iters = MatrixFreePcg::new(100, 1e-12)
            .solve(&op, &precond, &b, &mut x)
            .expect("PCG should converge on a small SPD system");
        assert!(iters <= 3, "CG on 3x3 SPD must converge in <= 3 iters");

        // Verify residual ‖b - A·x‖ is tiny.
        let mut ax = [0.0; 3];
        op.apply(&x, &mut ax);
        let res: f64 = b
            .iter()
            .zip(ax.iter())
            .map(|(bi, axi)| (bi - axi) * (bi - axi))
            .sum::<f64>()
            .sqrt();
        assert!(res < 1e-8, "residual too large: {res:.3e}");

        // Cross-check against the analytic solution (exact for this small system).
        // Solve by hand via the inverse: x = A^{-1} b.
        // det = 4*(4*3-1) - 1*(1*3-0) = 4*11 - 3 = 41.
        let det = 41.0;
        let x_exact = [
            (11.0 * b[0] - 3.0 * b[1] + 1.0 * b[2]) / det,
            (-3.0 * b[0] + 12.0 * b[1] - 4.0 * b[2]) / det,
            (1.0 * b[0] - 4.0 * b[1] + 15.0 * b[2]) / det,
        ];
        for i in 0..3 {
            assert!(
                (x[i] - x_exact[i]).abs() < 1e-10,
                "x[{i}] = {} expected {}",
                x[i],
                x_exact[i]
            );
        }
    }

    #[test]
    fn pcg_zero_rhs_returns_zero_iters() {
        let op = DenseOp {
            a: vec![4.0, 1.0, 0.0, 1.0, 4.0, 1.0, 0.0, 1.0, 3.0],
            n: 3,
        };
        let precond = IdentityPrecond { n: 3 };
        let b = [0.0; 3];
        let mut x = [0.0; 3];
        let iters = MatrixFreePcg::new(100, 1e-12)
            .solve(&op, &precond, &b, &mut x)
            .expect("zero rhs solves trivially");
        assert_eq!(iters, 0);
    }

    #[test]
    fn pcg_dimension_mismatch() {
        let op = DenseOp {
            a: vec![1.0],
            n: 1,
        };
        let precond = IdentityPrecond { n: 1 };
        let b = [1.0, 2.0];
        let mut x = [0.0];
        match MatrixFreePcg::new(10, 1e-10).solve(&op, &precond, &b, &mut x) {
            Err(MatrixFreeError::DimensionMismatch(1, 2)) => {}
            other => panic!("expected DimensionMismatch(1, 2), got {other:?}"),
        }
    }
}
