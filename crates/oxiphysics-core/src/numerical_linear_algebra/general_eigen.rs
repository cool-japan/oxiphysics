// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Eigenvalues of a general (non-symmetric) real square matrix.
//!
//! The matrix is reduced to real Schur form by nalgebra's Hessenberg +
//! implicit double-shift QR algorithm.  Unlike `DMatrix::complex_eigenvalues`,
//! which unwraps the decomposition and lets the QR iteration run without a
//! sweep limit, [`real_matrix_eigenvalues`] checks its input for non-finite
//! entries, bounds the number of sweeps and reports failure as `None`.
//! Callers (stability tests, characteristic-polynomial root finding) can
//! therefore never panic or hang on a diverged Jacobian.
//!
//! The double-shift iteration has no exceptional shifts, so it stalls on
//! matrices whose trailing blocks are invariant under it, the cyclic
//! permutation matrices that are the companion matrices of `x^n - c` being the
//! standard example.  When the first attempt exhausts its sweep budget, the
//! decomposition is retried on an orthogonally similar copy of the matrix
//! (same spectrum up to rounding, generic structure), which converges.

use nalgebra::{DMatrix, Schur};

/// Rotation angle (radians) of the first scrambling rotation.
const SCRAMBLE_ANGLE_OFFSET: f64 = 0.7;

/// Angle increment between consecutive scrambling rotations (about 74.5 degrees,
/// far from the multiples of pi/2 that would merely permute coordinates).
const SCRAMBLE_ANGLE_STEP: f64 = 1.3;

/// QR sweeps granted per eigenvalue before the decomposition is abandoned.
///
/// Double-shift QR needs two to three sweeps per eigenvalue in practice; a
/// budget of 100 only runs out on non-finite or pathologically scaled input.
const MAX_SWEEPS_PER_EIGENVALUE: usize = 100;

/// Eigenvalues `(re, im)` of the row-major `n x n` matrix `a`, sorted by real
/// part and then by imaginary part.
///
/// Real eigenvalues have `im == 0.0`; a complex-conjugate pair `re +/- i im`
/// appears as `(re, -im)` followed by `(re, +im)`.  `n == 0` yields an empty
/// spectrum.
///
/// Returns `None` when `a.len() != n * n`, when any entry is NaN or infinite,
/// or when the QR iteration does not converge within its sweep budget.
pub(crate) fn real_matrix_eigenvalues(a: &[f64], n: usize) -> Option<Vec<(f64, f64)>> {
    if n == 0 {
        return Some(Vec::new());
    }
    if a.len() != n.checked_mul(n)? || a.iter().any(|v| !v.is_finite()) {
        return None;
    }
    let matrix = DMatrix::from_row_slice(n, n, a);
    let budget = MAX_SWEEPS_PER_EIGENVALUE.checked_mul(n)?;
    let schur = match Schur::try_new(matrix.clone(), f64::EPSILON, budget) {
        Some(schur) => schur,
        None => Schur::try_new(scramble(matrix), f64::EPSILON, budget)?,
    };
    let mut spectrum: Vec<(f64, f64)> = schur
        .complex_eigenvalues()
        .iter()
        .map(|lambda| (lambda.re, lambda.im))
        .collect();
    if spectrum
        .iter()
        .any(|(re, im)| !re.is_finite() || !im.is_finite())
    {
        return None;
    }
    spectrum.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.total_cmp(&y.1)));
    Some(spectrum)
}

/// `Q A Qᵀ` for a fixed orthogonal `Q` made of Givens rotations of the
/// coordinate pairs `(i, i + 1)` and then `(i, i + 2)`.
///
/// The similarity preserves the spectrum up to rounding (orthogonal transforms
/// are backward stable) while destroying exact structure such as a cyclic
/// permutation pattern.
fn scramble(mut matrix: DMatrix<f64>) -> DMatrix<f64> {
    let n = matrix.nrows();
    let mut angle = SCRAMBLE_ANGLE_OFFSET;
    for stride in [1, 2] {
        for p in 0..n.saturating_sub(stride) {
            rotate_similarity(&mut matrix, p, p + stride, angle);
            angle += SCRAMBLE_ANGLE_STEP;
        }
    }
    matrix
}

/// `A <- G A Gᵀ` for the Givens rotation `G` of the coordinates `p < q` by `angle`.
fn rotate_similarity(matrix: &mut DMatrix<f64>, p: usize, q: usize, angle: f64) {
    let (sin, cos) = angle.sin_cos();
    for column in 0..matrix.ncols() {
        let (a, b) = (matrix[(p, column)], matrix[(q, column)]);
        matrix[(p, column)] = cos * a - sin * b;
        matrix[(q, column)] = sin * a + cos * b;
    }
    for row in 0..matrix.nrows() {
        let (a, b) = (matrix[(row, p)], matrix[(row, q)]);
        matrix[(row, p)] = cos * a - sin * b;
        matrix[(row, q)] = sin * a + cos * b;
    }
}

/// [`real_matrix_eigenvalues`] for a matrix stored as a slice of rows.
///
/// Returns `None` if the rows do not form a square matrix, in addition to the
/// failure modes of [`real_matrix_eigenvalues`].
pub(crate) fn real_matrix_eigenvalues_rows(rows: &[Vec<f64>]) -> Option<Vec<(f64, f64)>> {
    let n = rows.len();
    if rows.iter().any(|row| row.len() != n) {
        return None;
    }
    let flat: Vec<f64> = rows.iter().flat_map(|row| row.iter().copied()).collect();
    real_matrix_eigenvalues(&flat, n)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_spectrum(got: &[(f64, f64)], want: &[(f64, f64)], tol: f64) {
        assert_eq!(got.len(), want.len(), "spectrum {got:?} vs {want:?}");
        for (g, w) in got.iter().zip(want) {
            assert!(
                (g.0 - w.0).abs() <= tol && (g.1 - w.1).abs() <= tol,
                "spectrum {got:?} vs {want:?}"
            );
        }
    }

    #[test]
    fn empty_matrix_has_empty_spectrum() {
        assert_eq!(real_matrix_eigenvalues(&[], 0), Some(Vec::new()));
        assert_eq!(real_matrix_eigenvalues_rows(&[]), Some(Vec::new()));
    }

    #[test]
    fn diagonal_and_triangular_spectra_are_the_diagonal() {
        let a = [3.0, 1.0, 4.0, 0.0, -2.0, 5.0, 0.0, 0.0, 0.5];
        let eig = real_matrix_eigenvalues(&a, 3).expect("finite matrix");
        assert_spectrum(&eig, &[(-2.0, 0.0), (0.5, 0.0), (3.0, 0.0)], 1e-12);
    }

    #[test]
    fn rotation_block_gives_a_conjugate_pair() {
        // Eigenvalues of [[cos, -sin], [sin, cos]] are cos +/- i sin.
        let (c, s) = (0.6_f64, 0.8_f64);
        let eig = real_matrix_eigenvalues(&[c, -s, s, c], 2).expect("finite matrix");
        assert_spectrum(&eig, &[(c, -s), (c, s)], 1e-12);
    }

    #[test]
    fn companion_matrix_of_a_cubic_recovers_all_roots() {
        // (x - 1)(x^2 + 1) = x^3 - x^2 + x - 1: companion rows [0,1,0],[0,0,1],[1,-1,1].
        let rows = vec![
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
            vec![1.0, -1.0, 1.0],
        ];
        let eig = real_matrix_eigenvalues_rows(&rows).expect("finite matrix");
        assert_spectrum(&eig, &[(0.0, -1.0), (0.0, 1.0), (1.0, 0.0)], 1e-12);
    }

    #[test]
    fn symmetric_dominant_matrix_has_one_positive_eigenvalue() {
        // [[1,5,5],[5,1,5],[5,5,1]] has the spectrum {11, -4, -4}.
        let a = [1.0, 5.0, 5.0, 5.0, 1.0, 5.0, 5.0, 5.0, 1.0];
        let eig = real_matrix_eigenvalues(&a, 3).expect("finite matrix");
        assert_spectrum(&eig, &[(-4.0, 0.0), (-4.0, 0.0), (11.0, 0.0)], 1e-10);
    }

    #[test]
    fn invalid_input_is_reported_instead_of_panicking_or_hanging() {
        // Wrong length.
        assert!(real_matrix_eigenvalues(&[1.0, 2.0, 3.0], 2).is_none());
        // Non-finite entries.
        assert!(real_matrix_eigenvalues(&[1.0, f64::NAN, 0.0, 1.0], 2).is_none());
        assert!(real_matrix_eigenvalues(&[f64::INFINITY; 9], 3).is_none());
        // Ragged rows.
        let ragged = vec![vec![1.0, 2.0], vec![3.0]];
        assert!(real_matrix_eigenvalues_rows(&ragged).is_none());
        // Non-square row count.
        let wide = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]];
        assert!(real_matrix_eigenvalues_rows(&wide).is_none());
    }

    #[test]
    fn cyclic_permutation_matrices_converge() {
        // The companion matrix of x^n - 1 is the cyclic shift; its eigenvalues are
        // the n-th roots of unity.  The plain double-shift QR iteration leaves it
        // unchanged forever.
        for n in [3_usize, 5, 8, 12] {
            let mut a = vec![0.0; n * n];
            a[n - 1] = 1.0;
            for i in 1..n {
                a[i * n + (i - 1)] = 1.0;
            }
            let eig = real_matrix_eigenvalues(&a, n).expect("converges after scrambling");
            assert_eq!(eig.len(), n);
            let mut want: Vec<(f64, f64)> = (0..n)
                .map(|k| {
                    let angle = 2.0 * std::f64::consts::PI * k as f64 / n as f64;
                    (angle.cos(), angle.sin())
                })
                .collect();
            want.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.total_cmp(&y.1)));
            assert_spectrum(&eig, &want, 1e-12);
        }
    }

    #[test]
    fn scrambling_is_a_similarity_transform() {
        // Q A Qᵀ has the trace, the Frobenius norm and the determinant of A.
        let a = DMatrix::from_row_slice(
            4,
            4,
            &[
                1.0, 2.0, 0.0, -1.0, 0.5, -3.0, 4.0, 0.0, 2.0, 0.0, 1.0, 1.0, -1.0, 1.0, 0.0, 2.0,
            ],
        );
        let b = scramble(a.clone());
        assert!((a.trace() - b.trace()).abs() < 1e-12);
        assert!((a.norm() - b.norm()).abs() < 1e-12);
        assert!((a.determinant() - b.determinant()).abs() < 1e-10);
        assert!((&a - &b).norm() > 1e-3, "the matrix must actually change");
    }

    #[test]
    fn zero_matrix_has_a_zero_spectrum() {
        let eig = real_matrix_eigenvalues(&[0.0; 9], 3).expect("finite matrix");
        assert_spectrum(&eig, &[(0.0, 0.0), (0.0, 0.0), (0.0, 0.0)], 0.0);
    }
}
