// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Jacobi-preconditioned conjugate gradients for the sparse symmetric positive
//! definite systems of the parameterization algorithms.
//!
//! The solver stops on the *relative* residual `‖b − A x‖ / ‖b‖`, so it behaves
//! the same for a mesh measured in millimetres and one measured in kilometres.
//! (The generic `oxiphysics_core::sparse::pcg_solve` also gives up when
//! `pᵀ A p` drops below an absolute `f64::EPSILON`, which for a small right-hand
//! side stops the iteration long before the requested accuracy is reached.)

use crate::error::{Error, Result};
use oxiphysics_core::sparse::{CsrMatrix, diagonal_preconditioner};

/// Relative residual the iteration runs down to.
const TARGET_RESIDUAL: f64 = 1.0e-12;

/// Relative true residual above which a solve is reported as failed.
const ACCEPTABLE_RESIDUAL: f64 = 1.0e-9;

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

fn norm(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}

/// Solve `matrix · x = rhs` for a symmetric positive (semi-)definite `matrix`.
///
/// The iteration starts from zero.  A consistent *singular* system (a right-hand
/// side in the range of a positive semi-definite matrix) is solved as well, but
/// the preconditioner decides which of its solutions comes out: every iterate
/// satisfies `zᵀ D x = 0` for each null vector `z` of the matrix, where `D` is the
/// diagonal of `matrix`.  That is the minimal-norm solution in the
/// diagonal-weighted norm, not in the Euclidean one.  `operation` names the caller
/// in the error.
///
/// # Errors
/// * [`Error::DimensionMismatch`] if `matrix` is not square of the size of `rhs`;
/// * [`Error::DegenerateGeometry`] for a non-finite right-hand side;
/// * [`Error::ConvergenceFailure`] if the true relative residual of the final
///   iterate exceeds `1e-9` (an indefinite or very badly conditioned matrix).
pub(super) fn solve_spd(
    matrix: &CsrMatrix,
    rhs: &[f64],
    operation: &'static str,
) -> Result<Vec<f64>> {
    let n = rhs.len();
    if matrix.nrows != n || matrix.ncols != n {
        return Err(Error::dimension_mismatch(
            "matrix rows",
            matrix.nrows,
            "right-hand side",
            n,
        ));
    }
    let rhs_norm = norm(rhs);
    if !rhs_norm.is_finite() {
        return Err(Error::degenerate_geometry(
            "the right-hand side of a linear system is not finite",
        ));
    }
    if rhs_norm == 0.0 {
        return Ok(vec![0.0; n]);
    }

    let inverse_diagonal = diagonal_preconditioner(matrix);
    let max_iterations = 20 * n + 200;
    let mut x = vec![0.0; n];
    let mut residual = rhs.to_vec();
    let mut preconditioned: Vec<f64> = residual
        .iter()
        .zip(inverse_diagonal.iter())
        .map(|(r, m)| r * m)
        .collect();
    let mut direction = preconditioned.clone();
    let mut r_dot_z = dot(&residual, &preconditioned);
    let mut iterations = 0usize;
    while iterations < max_iterations && norm(&residual) > TARGET_RESIDUAL * rhs_norm {
        let image = matrix.matvec(&direction);
        let curvature = dot(&direction, &image);
        if curvature.is_nan() || curvature <= 0.0 {
            break;
        }
        let step = r_dot_z / curvature;
        for ((xi, ri), (di, ai)) in x
            .iter_mut()
            .zip(residual.iter_mut())
            .zip(direction.iter().zip(image.iter()))
        {
            *xi += step * di;
            *ri -= step * ai;
        }
        for ((zi, ri), mi) in preconditioned
            .iter_mut()
            .zip(residual.iter())
            .zip(inverse_diagonal.iter())
        {
            *zi = ri * mi;
        }
        let next_r_dot_z = dot(&residual, &preconditioned);
        let beta = next_r_dot_z / r_dot_z;
        for (di, zi) in direction.iter_mut().zip(preconditioned.iter()) {
            *di = zi + beta * *di;
        }
        r_dot_z = next_r_dot_z;
        iterations += 1;
    }

    // The recurrence residual drifts in finite precision: judge the true one.
    let image = matrix.matvec(&x);
    let true_residual = rhs
        .iter()
        .zip(image.iter())
        .map(|(b, ax)| (b - ax) * (b - ax))
        .sum::<f64>()
        .sqrt();
    if true_residual.is_nan() || true_residual > ACCEPTABLE_RESIDUAL * rhs_norm {
        return Err(Error::convergence_failure(
            operation,
            iterations,
            true_residual,
        ));
    }
    Ok(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tridiagonal `[-1, 2, -1]` matrix of size `n` (1-D Dirichlet Laplacian).
    fn laplacian_1d(n: usize) -> CsrMatrix {
        let mut row_ptr = vec![0];
        let mut col_idx = Vec::new();
        let mut values = Vec::new();
        for i in 0..n {
            if i > 0 {
                col_idx.push(i - 1);
                values.push(-1.0);
            }
            col_idx.push(i);
            values.push(2.0);
            if i + 1 < n {
                col_idx.push(i + 1);
                values.push(-1.0);
            }
            row_ptr.push(col_idx.len());
        }
        CsrMatrix {
            nrows: n,
            ncols: n,
            row_ptr,
            col_idx,
            values,
        }
    }

    #[test]
    fn solves_the_one_dimensional_laplacian_to_machine_precision() {
        let n = 200;
        let matrix = laplacian_1d(n);
        let rhs: Vec<f64> = (0..n).map(|i| ((i * 37) % 11) as f64 - 5.0).collect();
        let x = solve_spd(&matrix, &rhs, "test").unwrap();
        let image = matrix.matvec(&x);
        let error = rhs
            .iter()
            .zip(image.iter())
            .map(|(b, ax)| (b - ax).abs())
            .fold(0.0, f64::max);
        assert!(error < 1e-9, "max residual {error}");
    }

    #[test]
    fn accuracy_is_independent_of_the_scale_of_the_right_hand_side() {
        // `pcg_solve` of oxiphysics-core stops early for tiny right-hand sides.
        let n = 60;
        let matrix = laplacian_1d(n);
        let unit: Vec<f64> = (0..n).map(|i| (i as f64 * 0.37).sin()).collect();
        let reference = solve_spd(&matrix, &unit, "test").unwrap();
        for scale in [1e-12, 1e-6, 1e6, 1e12] {
            let rhs: Vec<f64> = unit.iter().map(|b| b * scale).collect();
            let x = solve_spd(&matrix, &rhs, "test").unwrap();
            for (a, b) in x.iter().zip(reference.iter()) {
                assert!((a - b * scale).abs() <= 1e-9 * scale.abs() * (1.0 + b.abs()));
            }
        }
    }

    #[test]
    fn zero_right_hand_side_gives_the_zero_solution() {
        let matrix = laplacian_1d(5);
        assert_eq!(solve_spd(&matrix, &[0.0; 5], "test").unwrap(), vec![0.0; 5]);
        let empty = CsrMatrix::new(0, 0);
        assert!(solve_spd(&empty, &[], "test").unwrap().is_empty());
    }

    #[test]
    fn consistent_singular_systems_are_solved_with_the_kernel_part_fixed_by_the_diagonal() {
        // Pure Neumann Laplacian: the constants span the kernel; its diagonal is
        // (1, 2, 2, 2, 2, 1).
        let mut matrix = laplacian_1d(6);
        matrix.values[0] = 1.0;
        let last_diagonal = matrix.values.len() - 1;
        matrix.values[last_diagonal] = 1.0;
        let rhs = [1.0, -1.0, 0.0, 0.0, 0.5, -0.5];
        let x = solve_spd(&matrix, &rhs, "test").unwrap();
        let image = matrix.matvec(&x);
        for (b, ax) in rhs.iter().zip(image.iter()) {
            assert!((b - ax).abs() < 1e-9, "A x = {image:?}");
        }
        // The null vector z = (1, ..., 1) is orthogonal to x in the inner product
        // weighted by the diagonal of the matrix (not in the Euclidean one).
        let diagonal = [1.0, 2.0, 2.0, 2.0, 2.0, 1.0];
        let weighted: f64 = x.iter().zip(diagonal.iter()).map(|(xi, di)| xi * di).sum();
        assert!(weighted.abs() < 1e-9, "z^T D x = {weighted}");
        // Adding any multiple of the null vector gives another solution: the
        // unweighted mean is not constrained.
        let shifted: Vec<f64> = x.iter().map(|xi| xi + 0.75).collect();
        let shifted_image = matrix.matvec(&shifted);
        for (b, ax) in rhs.iter().zip(shifted_image.iter()) {
            assert!((b - ax).abs() < 1e-9);
        }
    }

    #[test]
    fn invalid_systems_are_reported() {
        let matrix = laplacian_1d(4);
        assert!(matches!(
            solve_spd(&matrix, &[1.0, 2.0, 3.0], "test"),
            Err(Error::DimensionMismatch { .. })
        ));
        assert!(matches!(
            solve_spd(&matrix, &[1.0, f64::NAN, 0.0, 0.0], "test"),
            Err(Error::DegenerateGeometry { .. })
        ));
        // Indefinite matrix: diag(1, -1).
        let indefinite = CsrMatrix {
            nrows: 2,
            ncols: 2,
            row_ptr: vec![0, 1, 2],
            col_idx: vec![0, 1],
            values: vec![1.0, -1.0],
        };
        let failure = solve_spd(&indefinite, &[1.0, 1.0], "indefinite test");
        assert!(matches!(
            failure,
            Err(Error::ConvergenceFailure {
                operation: "indefinite test",
                ..
            })
        ));
    }
}
