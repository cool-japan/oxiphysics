// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Dense LU factorisation with partial pivoting and explicit singularity detection.
//!
//! The implicit time integrators of [`crate::dynamics`] and the dense eigen-solvers of
//! [`crate::eigen`] solve small dense systems whose matrix may legitimately be (numerically)
//! singular, for example when a shift sits exactly on an eigenvalue or when a structure
//! carries rigid-body modes.  Those callers must be *told* about the singularity instead of
//! receiving a silently perturbed or zeroed solution, which is what this factorisation
//! guarantees: [`LuFactor::new`] returns `None` as soon as a pivot is not larger than the
//! rounding level of the matrix.

/// LU factors `P·A = L·U` of a square matrix, computed with partial (row) pivoting.
#[derive(Debug, Clone)]
pub(crate) struct LuFactor {
    /// Matrix dimension.
    n: usize,
    /// Row-major packed factors: the strict lower triangle holds `L` (its unit diagonal is
    /// implicit), the remaining entries hold `U`.
    lu: Vec<f64>,
    /// `perm[i]` is the index of the original row that was moved to row `i`.
    perm: Vec<usize>,
}

impl LuFactor {
    /// Factor the square matrix `a`.
    ///
    /// Returns `None` when the matrix is not square, contains a non-finite entry, or is
    /// singular to working precision, i.e. when a pivot is not larger than
    /// `n · ε · max|a_ij|`.  An empty matrix factors trivially.
    pub(crate) fn new(a: &[Vec<f64>]) -> Option<Self> {
        let n = a.len();
        let mut lu: Vec<f64> = Vec::with_capacity(n * n);
        let mut max_abs = 0.0_f64;
        for row in a {
            if row.len() != n {
                return None;
            }
            for &x in row {
                if !x.is_finite() {
                    return None;
                }
                max_abs = max_abs.max(x.abs());
                lu.push(x);
            }
        }
        let tol = n as f64 * f64::EPSILON * max_abs;
        let mut perm: Vec<usize> = (0..n).collect();
        for k in 0..n {
            let mut pivot_row = k;
            let mut best = lu[k * n + k].abs();
            for r in (k + 1)..n {
                let candidate = lu[r * n + k].abs();
                if candidate > best {
                    best = candidate;
                    pivot_row = r;
                }
            }
            if best.is_nan() || best <= tol {
                return None;
            }
            if pivot_row != k {
                for j in 0..n {
                    lu.swap(k * n + j, pivot_row * n + j);
                }
                perm.swap(k, pivot_row);
            }
            let pivot = lu[k * n + k];
            for r in (k + 1)..n {
                let factor = lu[r * n + k] / pivot;
                lu[r * n + k] = factor;
                if factor != 0.0 {
                    for j in (k + 1)..n {
                        let update = factor * lu[k * n + j];
                        lu[r * n + j] -= update;
                    }
                }
            }
        }
        Some(Self { n, lu, perm })
    }

    /// Solve `A·x = b` with the stored factors.
    ///
    /// `b` must have the dimension of the factored matrix; missing entries are read as zero.
    pub(crate) fn solve(&self, b: &[f64]) -> Vec<f64> {
        let n = self.n;
        let mut x: Vec<f64> = self
            .perm
            .iter()
            .map(|&row| b.get(row).copied().unwrap_or(0.0))
            .collect();
        // Forward substitution with the unit lower triangle: L·y = P·b.
        for i in 1..n {
            let (solved, rest) = x.split_at_mut(i);
            let l_row = &self.lu[i * n..i * n + i];
            let dot: f64 = l_row.iter().zip(solved.iter()).map(|(l, y)| l * y).sum();
            rest[0] -= dot;
        }
        // Back substitution with the upper triangle: U·x = y.
        for i in (0..n).rev() {
            let (head, solved) = x.split_at_mut(i + 1);
            let u_row = &self.lu[i * n + i + 1..(i + 1) * n];
            let dot: f64 = u_row.iter().zip(solved.iter()).map(|(u, xj)| u * xj).sum();
            head[i] = (head[i] - dot) / self.lu[i * n + i];
        }
        x
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn matvec(a: &[Vec<f64>], x: &[f64]) -> Vec<f64> {
        a.iter()
            .map(|row| row.iter().zip(x.iter()).map(|(aij, xj)| aij * xj).sum())
            .collect()
    }

    #[test]
    fn solves_a_system_that_needs_row_pivoting() {
        // The leading entry is zero, so an unpivoted elimination would divide by zero.
        let a = vec![
            vec![0.0, 2.0, 1.0],
            vec![1.0, 1.0, 0.0],
            vec![3.0, 0.0, 4.0],
        ];
        let x_true = vec![1.0, -2.0, 0.5];
        let b = matvec(&a, &x_true);
        let lu = LuFactor::new(&a).expect("regular matrix");
        let x = lu.solve(&b);
        for (xi, ti) in x.iter().zip(x_true.iter()) {
            assert!((xi - ti).abs() < 1e-13, "x = {x:?}");
        }
    }

    #[test]
    fn solves_nonsymmetric_and_indefinite_systems() {
        let a = vec![
            vec![2.0, -3.0, 1.0, 0.0],
            vec![-1.0, -4.0, 2.0, 5.0],
            vec![0.5, 1.0, -6.0, 2.0],
            vec![4.0, 0.0, 1.0, -2.0],
        ];
        let x_true = vec![0.3, -1.2, 2.5, 0.7];
        let b = matvec(&a, &x_true);
        let lu = LuFactor::new(&a).expect("regular matrix");
        let x = lu.solve(&b);
        for (xi, ti) in x.iter().zip(x_true.iter()) {
            assert!((xi - ti).abs() < 1e-12, "x = {x:?}");
        }
    }

    #[test]
    fn reports_singular_matrices() {
        let rank_deficient = vec![vec![1.0, 2.0], vec![2.0, 4.0]];
        assert!(LuFactor::new(&rank_deficient).is_none());
        let zero = vec![vec![0.0; 3]; 3];
        assert!(LuFactor::new(&zero).is_none());
    }

    #[test]
    fn singularity_test_is_relative_to_the_matrix_scale() {
        // The same rank-one matrix at a huge scale must still be recognised as singular,
        // and a well-conditioned matrix at a tiny scale must still factor.
        let big = vec![vec![1.0e12, 2.0e12], vec![2.0e12, 4.0e12]];
        assert!(LuFactor::new(&big).is_none());
        let tiny = vec![vec![2.0e-12, 1.0e-12], vec![1.0e-12, 3.0e-12]];
        let lu = LuFactor::new(&tiny).expect("regular matrix at a small scale");
        let x = lu.solve(&[5.0e-12, 7.0e-12]);
        assert!(
            (x[0] - 1.6).abs() < 1e-12 && (x[1] - 1.8).abs() < 1e-12,
            "{x:?}"
        );
    }

    #[test]
    fn rejects_non_square_and_non_finite_input() {
        assert!(LuFactor::new(&[vec![1.0, 2.0]]).is_none());
        assert!(LuFactor::new(&[vec![1.0, f64::NAN], vec![0.0, 1.0]]).is_none());
        assert!(LuFactor::new(&[vec![f64::INFINITY]]).is_none());
    }

    #[test]
    fn empty_matrix_factors_trivially() {
        let lu = LuFactor::new(&[]).expect("empty matrix");
        assert!(lu.solve(&[]).is_empty());
    }
}
