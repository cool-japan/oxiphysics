// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! θ-method time integration of the semi-discrete heat equation `C Ṫ + K T = f`.
//!
//! The solver advances nodal temperatures with explicit Euler (θ = 0),
//! Crank–Nicolson (θ = ½) or implicit Euler (θ = 1).  Implicit variants solve the
//! full sparse system, off-diagonal conduction included.

use crate::analysis::FemError;
use crate::solvers::{bicgstab_sparse, pcg_solve};
use crate::sparse::CsrMatrix;

/// Transient solver using the generalised θ-method.
///
/// θ = 0  → explicit Euler.
/// θ = 0.5 → Crank–Nicolson (second-order accurate, unconditionally stable).
/// θ = 1  → implicit Euler.
///
/// Every step advances the semi-discrete system `C Ṫ + K T = f` (lumped `C`,
/// constant load `f`) by
///
/// ```text
/// (C + θ dt K) T_{n+1} = (C − (1−θ) dt K) T_n + dt f.
/// ```
///
/// For θ > 0 the *full* sparse system — off-diagonal conduction included — is
/// solved with Jacobi-preconditioned conjugate gradients (BiCGSTAB when `K` is
/// not symmetric), so a field with `K T = f` is a fixed point of the scheme and
/// a uniform temperature is preserved by an adiabatic mesh.
pub struct TransientThermalSolver {
    /// θ parameter in \[0, 1\].
    pub theta: f64,
    /// Time step \[s\].
    pub dt: f64,
}

/// Time-invariant data of a θ-method step: lumped capacitance and (for θ > 0)
/// the sparse left-hand side `C + θ dt K` with its Jacobi preconditioner.
struct ThetaSystem {
    /// Lumped (row-sum) capacitance of every node.
    c_lumped: Vec<f64>,
    /// `C + θ dt K` in CSR form; `None` for the explicit scheme (θ = 0).
    lhs: Option<CsrMatrix>,
    /// Inverse diagonal of `lhs` (Jacobi preconditioner).
    jacobi: Vec<f64>,
    /// Whether `lhs` is symmetric (conjugate gradients) or not (BiCGSTAB).
    symmetric: bool,
}

impl TransientThermalSolver {
    /// Create a solver with explicit Euler (θ = 0).
    pub fn explicit_euler(dt: f64) -> Self {
        Self { theta: 0.0_f64, dt }
    }

    /// Create a Crank–Nicolson solver (θ = 0.5).
    pub fn crank_nicolson(dt: f64) -> Self {
        Self { theta: 0.5_f64, dt }
    }

    /// Create an implicit Euler solver (θ = 1).
    pub fn implicit_euler(dt: f64) -> Self {
        Self { theta: 1.0_f64, dt }
    }

    /// Advance one time step.
    ///
    /// Solves `(C + θ dt K) T_{n+1} = (C − (1−θ) dt K) T_n + dt f`, with the
    /// lumped capacitance `C` taken from the row sums of `m`.  Internally the
    /// equivalent increment form `(C + θ dt K) ΔT = dt (f − K T_n)` is solved, so
    /// a state with `K T_n = f` is reproduced exactly and no accuracy is lost to
    /// large absolute temperatures.
    ///
    /// # Errors
    /// * [`FemError::InvalidInput`] — `dt` is not positive and finite, θ is
    ///   outside \[0, 1\], `k`/`m` are not square `n × n` matrices matching
    ///   `temps` and `f`, or an input contains a non-finite value.
    /// * [`FemError::SingularSystem`] — `C + θ dt K` has a non-positive diagonal.
    /// * [`FemError::SolverNotConverged`] — the iterative solve of the implicit
    ///   system did not reach a relative residual of `1e-12`.
    pub fn step(
        &self,
        temps: &[f64],
        k: &[Vec<f64>],
        m: &[Vec<f64>],
        f: &[f64],
    ) -> Result<Vec<f64>, FemError> {
        let system = self.prepare(temps.len(), k, m)?;
        self.advance(&system, temps, k, f)
    }

    /// Run `n_steps` time integration steps, returning temperature history.
    ///
    /// The history has `n_steps + 1` entries (the initial state first).  The
    /// system matrix is assembled once and reused by every step.
    ///
    /// # Errors
    /// Same conditions as [`step`](Self::step).
    pub fn integrate(
        &self,
        t0: &[f64],
        k: &[Vec<f64>],
        m: &[Vec<f64>],
        f: &[f64],
        n_steps: usize,
    ) -> Result<Vec<Vec<f64>>, FemError> {
        let system = self.prepare(t0.len(), k, m)?;
        let mut history = Vec::with_capacity(n_steps + 1);
        let mut temps = t0.to_vec();
        history.push(temps.clone());
        for _ in 0..n_steps {
            temps = self.advance(&system, &temps, k, f)?;
            history.push(temps.clone());
        }
        Ok(history)
    }

    /// Validate the parameters and the matrices and build the step data.
    fn prepare(&self, n: usize, k: &[Vec<f64>], m: &[Vec<f64>]) -> Result<ThetaSystem, FemError> {
        let dt = self.dt;
        let theta = self.theta;
        if !dt.is_finite() || dt <= 0.0 {
            return Err(FemError::InvalidInput(format!(
                "time step must be positive and finite, got {dt}"
            )));
        }
        if !(0.0..=1.0).contains(&theta) {
            return Err(FemError::InvalidInput(format!(
                "theta must lie in [0, 1], got {theta}"
            )));
        }
        let is_square = |a: &[Vec<f64>]| a.len() == n && a.iter().all(|row| row.len() == n);
        if !is_square(k) || !is_square(m) {
            return Err(FemError::InvalidInput(format!(
                "K and M must be square {n} x {n} matrices matching the temperature vector"
            )));
        }
        if k.iter()
            .chain(m.iter())
            .any(|row| row.iter().any(|v| !v.is_finite()))
        {
            return Err(FemError::InvalidInput(
                "K and M must contain only finite values".to_string(),
            ));
        }

        // Lumped capacitance: sum of each row of M
        let c_lumped: Vec<f64> = m
            .iter()
            .map(|row| row.iter().sum::<f64>().max(1e-30))
            .collect();

        if theta == 0.0_f64 {
            return Ok(ThetaSystem {
                c_lumped,
                lhs: None,
                jacobi: Vec::new(),
                symmetric: true,
            });
        }

        // Left-hand side C + theta*dt*K in CSR form (the diagonal is always stored).
        let mut row_ptr = Vec::with_capacity(n + 1);
        let mut col_indices = Vec::new();
        let mut values = Vec::new();
        let mut jacobi = Vec::with_capacity(n);
        row_ptr.push(0);
        for (i, (k_row, &c_i)) in k.iter().zip(c_lumped.iter()).enumerate() {
            for (j, &k_ij) in k_row.iter().enumerate() {
                let mut a_ij = theta * dt * k_ij;
                if i == j {
                    a_ij += c_i;
                    if a_ij <= 0.0 {
                        return Err(FemError::SingularSystem);
                    }
                    jacobi.push(1.0 / a_ij);
                }
                if i == j || a_ij != 0.0 {
                    col_indices.push(j);
                    values.push(a_ij);
                }
            }
            row_ptr.push(col_indices.len());
        }
        let symmetric = (0..n).all(|i| {
            (0..i).all(|j| {
                let (a, b) = (k[i][j], k[j][i]);
                (a - b).abs() <= 1.0e-10 * a.abs().max(b.abs())
            })
        });
        Ok(ThetaSystem {
            c_lumped,
            lhs: Some(CsrMatrix {
                row_ptr,
                col_indices,
                values,
                nrows: n,
                ncols: n,
            }),
            jacobi,
            symmetric,
        })
    }

    /// Advance `temps` by one step with the prepared system.
    fn advance(
        &self,
        system: &ThetaSystem,
        temps: &[f64],
        k: &[Vec<f64>],
        f: &[f64],
    ) -> Result<Vec<f64>, FemError> {
        let n = system.c_lumped.len();
        if temps.len() != n || f.len() != n {
            return Err(FemError::InvalidInput(format!(
                "temperature ({}) and load ({}) vectors must have {n} entries",
                temps.len(),
                f.len()
            )));
        }
        if temps.iter().chain(f.iter()).any(|v| !v.is_finite()) {
            return Err(FemError::InvalidInput(
                "temperatures and loads must be finite".to_string(),
            ));
        }
        let dt = self.dt;

        // dt * (f - K T_n): the change of the right-hand side that drives the step.
        let drive: Vec<f64> = k
            .iter()
            .zip(f.iter())
            .map(|(k_row, &f_i)| {
                let kt_i: f64 = k_row
                    .iter()
                    .zip(temps.iter())
                    .map(|(k_ij, t_j)| k_ij * t_j)
                    .sum();
                dt * (f_i - kt_i)
            })
            .collect();

        let delta = match &system.lhs {
            // Explicit: C (T_{n+1} - T_n) = dt (f - K T_n)
            None => drive
                .iter()
                .zip(system.c_lumped.iter())
                .map(|(d, c)| d / c)
                .collect(),
            // Implicit: (C + theta*dt*K) (T_{n+1} - T_n) = dt (f - K T_n)
            Some(lhs) => Self::solve_increment(system, lhs, &drive)?,
        };
        Ok(temps.iter().zip(delta.iter()).map(|(t, d)| t + d).collect())
    }

    /// Solve `lhs · ΔT = drive` (scaled so the solver tolerances are relative).
    fn solve_increment(
        system: &ThetaSystem,
        lhs: &CsrMatrix,
        drive: &[f64],
    ) -> Result<Vec<f64>, FemError> {
        let scale = drive.iter().fold(0.0_f64, |acc, v| acc.max(v.abs()));
        if scale == 0.0_f64 {
            // Nothing drives the field: an exact fixed point of the scheme.
            return Ok(vec![0.0_f64; drive.len()]);
        }
        let rhs: Vec<f64> = drive.iter().map(|v| v / scale).collect();
        let max_iter = (10 * drive.len()).max(100);
        let (y, stats) = if system.symmetric {
            pcg_solve(lhs, &rhs, &system.jacobi, max_iter, 1.0e-12)
        } else {
            bicgstab_sparse(lhs, &rhs, max_iter, 1.0e-12)
        };
        if !stats.converged || y.iter().any(|v| !v.is_finite()) {
            return Err(FemError::SolverNotConverged);
        }
        Ok(y.into_iter().map(|v| v * scale).collect())
    }
}
