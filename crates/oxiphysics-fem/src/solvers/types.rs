//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
#[allow(unused_imports)]
use super::functions_2::*;

/// Conjugate Gradient solver for symmetric positive definite systems.
#[allow(dead_code)]
pub struct ConjugateGradient {
    /// Maximum number of iterations.
    pub max_iter: usize,
    /// Convergence tolerance on the residual norm.
    pub tol: f64,
}
#[allow(dead_code)]
impl ConjugateGradient {
    /// Create a new CG solver.
    pub fn new(max_iter: usize, tol: f64) -> Self {
        Self { max_iter, tol }
    }
    /// Solve the dense system `A x = b` where `A` is stored in row-major order.
    ///
    /// `n` is the side length.  Returns the solution vector.
    pub fn solve(&self, a: &[f64], b: &[f64], n: usize) -> Vec<f64> {
        assert_eq!(a.len(), n * n);
        assert_eq!(b.len(), n);
        let mut x = vec![0.0; n];
        let mut r = b.to_vec();
        let mut p = r.clone();
        let mut rs_old = dot_slice(&r, &r);
        if rs_old.sqrt() < self.tol {
            return x;
        }
        for _ in 0..self.max_iter {
            let ap = dense_matvec(a, &p, n);
            let p_ap = dot_slice(&p, &ap);
            if p_ap.abs() < 1e-60 {
                break;
            }
            let alpha = rs_old / p_ap;
            for i in 0..n {
                x[i] += alpha * p[i];
                r[i] -= alpha * ap[i];
            }
            let rs_new = dot_slice(&r, &r);
            if rs_new.sqrt() < self.tol {
                break;
            }
            let beta = rs_new / rs_old;
            for i in 0..n {
                p[i] = r[i] + beta * p[i];
            }
            rs_old = rs_new;
        }
        x
    }
}
/// Error type returned by dense solvers.
#[derive(Debug, Clone, PartialEq)]
pub enum SolverError {
    /// Matrix is not positive definite (negative pivot encountered).
    NotPositiveDefinite,
    /// Dimension mismatch.
    DimensionMismatch,
}
/// Dense Cholesky factorisation (lower-triangular `L` such that `A = L Lᵀ`).
#[allow(dead_code)]
pub struct CholeskyDense {
    /// Dimension of the system.
    pub n: usize,
}
#[allow(dead_code)]
impl CholeskyDense {
    /// Create a new Cholesky solver.
    pub fn new(n: usize) -> Self {
        Self { n }
    }
    /// Factorise `A` (row-major, `n × n`, SPD) into lower-triangular `L`.
    ///
    /// Returns `Ok(L_flat)` where `L_flat` is stored in row-major order.
    pub fn factorize(&self, a: &[f64]) -> Result<Vec<f64>, SolverError> {
        let n = self.n;
        if a.len() != n * n {
            return Err(SolverError::DimensionMismatch);
        }
        let mut l = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..=i {
                let mut s = a[i * n + j];
                for k in 0..j {
                    s -= l[i * n + k] * l[j * n + k];
                }
                if i == j {
                    if s <= 0.0 {
                        return Err(SolverError::NotPositiveDefinite);
                    }
                    l[i * n + j] = s.sqrt();
                } else {
                    l[i * n + j] = s / l[j * n + j];
                }
            }
        }
        Ok(l)
    }
    /// Solve `A x = b` given the Cholesky factor `l` and the RHS `b`.
    ///
    /// Performs forward and backward substitution.
    pub fn solve(&self, l: &[f64], b: &[f64]) -> Vec<f64> {
        let n = self.n;
        let mut y = vec![0.0; n];
        for i in 0..n {
            let mut s = b[i];
            for k in 0..i {
                s -= l[i * n + k] * y[k];
            }
            y[i] = s / l[i * n + i];
        }
        let mut x = vec![0.0; n];
        for i in (0..n).rev() {
            let mut s = y[i];
            for k in (i + 1)..n {
                s -= l[k * n + i] * x[k];
            }
            x[i] = s / l[i * n + i];
        }
        x
    }
}
/// Convergence monitor that records residual history during iterative solving.
///
/// Attach to any iterative solver to track convergence and detect stagnation
/// or divergence.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ConvergenceMonitor {
    /// Residual history (one entry per iteration).
    pub history: Vec<f64>,
    /// Convergence tolerance.
    pub tol: f64,
    /// Maximum allowed iterations.
    pub max_iter: usize,
}
#[allow(dead_code)]
impl ConvergenceMonitor {
    /// Create a new convergence monitor.
    pub fn new(tol: f64, max_iter: usize) -> Self {
        Self {
            history: Vec::new(),
            tol,
            max_iter,
        }
    }
    /// Record a residual value for the current iteration.
    pub fn record(&mut self, residual: f64) {
        self.history.push(residual);
    }
    /// Return `true` if the last recorded residual is below the tolerance.
    pub fn converged(&self) -> bool {
        self.history.last().is_some_and(|&r| r < self.tol)
    }
    /// Return `true` if the solver is stagnating (last two residuals differ
    /// by less than 1e-14 relative).
    pub fn stagnated(&self) -> bool {
        if self.history.len() < 2 {
            return false;
        }
        let n = self.history.len();
        let r0 = self.history[n - 2];
        let r1 = self.history[n - 1];
        let scale = r0.abs().max(1e-60);
        (r1 - r0).abs() / scale < 1e-14
    }
    /// Return the convergence rate (ratio of successive residuals), or `None`
    /// if fewer than 2 iterations have been recorded.
    pub fn convergence_rate(&self) -> Option<f64> {
        let n = self.history.len();
        if n < 2 || self.history[n - 2].abs() < 1e-60 {
            None
        } else {
            Some(self.history[n - 1] / self.history[n - 2])
        }
    }
    /// Return the number of iterations recorded.
    pub fn iterations(&self) -> usize {
        self.history.len()
    }
}
/// Gauss-Seidel (successive over-relaxation) solver for dense systems.
#[allow(dead_code)]
pub struct GaussSeidel {
    /// Maximum number of iterations.
    pub max_iter: usize,
    /// Convergence tolerance.
    pub tol: f64,
    /// Relaxation parameter (1.0 = standard Gauss-Seidel, >1 = over-relaxation).
    pub omega: f64,
}
#[allow(dead_code)]
impl GaussSeidel {
    /// Create a new Gauss-Seidel solver.
    pub fn new(max_iter: usize, tol: f64, omega: f64) -> Self {
        Self {
            max_iter,
            tol,
            omega,
        }
    }
    /// Solve the dense row-major system `A x = b` of size `n`.
    ///
    /// Returns the solution vector.
    pub fn solve(&self, a: &[f64], b: &[f64], n: usize) -> Vec<f64> {
        let mut x = vec![0.0; n];
        for _ in 0..self.max_iter {
            let mut max_change = 0.0_f64;
            for i in 0..n {
                let mut sigma = 0.0;
                for j in 0..n {
                    if j != i {
                        sigma += a[i * n + j] * x[j];
                    }
                }
                let diag = a[i * n + i];
                if diag.abs() < 1e-60 {
                    continue;
                }
                let x_new = (1.0 - self.omega) * x[i] + self.omega * (b[i] - sigma) / diag;
                max_change = max_change.max((x_new - x[i]).abs());
                x[i] = x_new;
            }
            if max_change < self.tol {
                break;
            }
        }
        x
    }
    /// Solve and return `(solution, stats)`.
    pub fn solve_with_stats(&self, a: &[f64], b: &[f64], n: usize) -> (Vec<f64>, SolverStats) {
        let mut x = vec![0.0; n];
        let mut iters = 0;
        let mut converged = false;
        for _it in 0..self.max_iter {
            iters += 1;
            let mut max_change = 0.0_f64;
            for i in 0..n {
                let mut sigma = 0.0;
                for j in 0..n {
                    if j != i {
                        sigma += a[i * n + j] * x[j];
                    }
                }
                let diag = a[i * n + i];
                if diag.abs() < 1e-60 {
                    continue;
                }
                let x_new = (1.0 - self.omega) * x[i] + self.omega * (b[i] - sigma) / diag;
                max_change = max_change.max((x_new - x[i]).abs());
                x[i] = x_new;
            }
            if max_change < self.tol {
                converged = true;
                break;
            }
        }
        let ax = dense_matvec(a, &x, n);
        let r: Vec<f64> = b.iter().zip(ax.iter()).map(|(bi, ai)| bi - ai).collect();
        let residual = dot_slice(&r, &r).sqrt();
        (
            x,
            SolverStats {
                iterations: iters,
                residual,
                converged,
            },
        )
    }
}
/// Statistics returned by iterative solvers.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct SolverStats {
    /// Number of iterations performed.
    pub iterations: usize,
    /// Final residual norm ‖r‖₂.
    pub residual: f64,
    /// True if the solver converged within the tolerance.
    pub converged: bool,
}
