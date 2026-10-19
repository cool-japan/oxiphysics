//! Regression tests for the Krylov eigensolvers (Arnoldi, Lanczos), LOBPCG,
//! TFQMR, the matrix logarithm and the tridiagonal / Hessenberg eigenvalue
//! routines.  Every expected value is either a closed form or an independent
//! residual check, never a number copied from the implementation under test.

use std::f64::consts::PI;

use super::{EigenSolver, IterativeSolver, MatrixFunctions, SparseCSR};
use crate::numerical_methods::symmetric_eigen_n;

fn laplacian(n: usize) -> SparseCSR {
    SparseCSR::laplacian_1d(n)
}

/// Exact eigenvalue `2 - 2 cos(j pi / (n + 1))` of the 1-D Laplacian.
fn laplacian_eigenvalue(j: usize, n: usize) -> f64 {
    2.0 - 2.0 * (j as f64 * PI / (n as f64 + 1.0)).cos()
}

fn scaled(a: &SparseCSR, factor: f64) -> SparseCSR {
    SparseCSR {
        nrows: a.nrows,
        ncols: a.ncols,
        row_ptr: a.row_ptr.clone(),
        col_idx: a.col_idx.clone(),
        values: a.values.iter().map(|v| v * factor).collect(),
    }
}

/// Sparse matrix from a dense row list (zeros are not stored).
fn from_dense(rows: &[Vec<f64>]) -> SparseCSR {
    let n = rows.len();
    let (mut r, mut c, mut v) = (Vec::new(), Vec::new(), Vec::new());
    for (i, row) in rows.iter().enumerate() {
        for (j, &value) in row.iter().enumerate() {
            if value != 0.0 {
                r.push(i);
                c.push(j);
                v.push(value);
            }
        }
    }
    SparseCSR::from_triplets(n, n, &r, &c, &v)
}

fn diagonal_matrix(entries: &[f64]) -> SparseCSR {
    let n = entries.len();
    let index: Vec<usize> = (0..n).collect();
    SparseCSR::from_triplets(n, n, &index, &index, entries)
}

/// Non-symmetric tridiagonal Toeplitz matrix with diagonal 2, sub-diagonal
/// `-1 - c` and super-diagonal `-1 + c`.  For `|c| < 1` its eigenvalues are real:
/// `2 + 2 sqrt(1 - c^2) cos(j pi / (n + 1))`.
fn convection_diffusion(n: usize, c: f64) -> SparseCSR {
    let (mut r, mut col, mut v) = (Vec::new(), Vec::new(), Vec::new());
    for i in 0..n {
        r.push(i);
        col.push(i);
        v.push(2.0);
        if i > 0 {
            r.push(i);
            col.push(i - 1);
            v.push(-1.0 - c);
        }
        if i + 1 < n {
            r.push(i);
            col.push(i + 1);
            v.push(-1.0 + c);
        }
    }
    SparseCSR::from_triplets(n, n, &r, &col, &v)
}

fn convection_diffusion_eigenvalue(j: usize, n: usize, c: f64) -> f64 {
    2.0 + 2.0 * (1.0 - c * c).sqrt() * (j as f64 * PI / (n as f64 + 1.0)).cos()
}

fn residual_norm(a: &SparseCSR, lambda: f64, v: &[f64]) -> f64 {
    a.matvec(v)
        .iter()
        .zip(v.iter())
        .map(|(av, x)| (av - lambda * x).powi(2))
        .sum::<f64>()
        .sqrt()
}

fn norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

fn inner(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

fn true_residual(a: &SparseCSR, x: &[f64], b: &[f64]) -> f64 {
    a.matvec(x)
        .iter()
        .zip(b.iter())
        .map(|(ax, bi)| (ax - bi).powi(2))
        .sum::<f64>()
        .sqrt()
}

// ---------------------------------------------------------------------------
// F0793: tridiagonal eigenvalues
// ---------------------------------------------------------------------------

#[test]
fn tridiag_eigenvalues_three_by_three_laplacian() {
    // The old fake Jacobi sweep returned [2, 2, 2].
    let ev = EigenSolver::tridiag_eigenvalues(&[2.0, 2.0, 2.0], &[-1.0, -1.0]);
    let s2 = std::f64::consts::SQRT_2;
    let expected = [2.0 - s2, 2.0, 2.0 + s2];
    assert_eq!(ev.len(), 3);
    for (got, want) in ev.iter().zip(expected.iter()) {
        assert!((got - want).abs() < 1e-12, "got {got}, want {want}");
    }
}

#[test]
fn tridiag_eigenvalues_varying_diagonal() {
    // numpy.linalg.eigvalsh of [[2,-1,0],[-1,3,.5],[0,.5,1]]; the old code
    // produced [-6.97, -1.62, 14.59].
    let ev = EigenSolver::tridiag_eigenvalues(&[2.0, 3.0, 1.0], &[-1.0, 0.5]);
    let expected = [0.8138593383654928, 1.5, 3.686140661634507];
    for (got, want) in ev.iter().zip(expected.iter()) {
        assert!((got - want).abs() < 1e-12, "got {got}, want {want}");
    }
    assert_eq!(ev.len(), 3);
}

#[test]
fn tridiag_eigenvalues_tolerates_surplus_and_missing_off_diagonals() {
    // A Lanczos run may hand over one off-diagonal too many: it must be ignored.
    let with_surplus = EigenSolver::tridiag_eigenvalues(&[2.0, 2.0, 2.0], &[-1.0, -1.0, 99.0]);
    let exact = EigenSolver::tridiag_eigenvalues(&[2.0, 2.0, 2.0], &[-1.0, -1.0]);
    assert_eq!(with_surplus, exact);
    // A missing entry decouples the trailing block: eigenvalues {1, 3} and {2}.
    let decoupled = EigenSolver::tridiag_eigenvalues(&[2.0, 2.0, 2.0], &[-1.0]);
    for (got, want) in decoupled.iter().zip([1.0, 2.0, 3.0].iter()) {
        assert!((got - want).abs() < 1e-12, "got {got}, want {want}");
    }
}

#[test]
fn tridiag_eigenvalues_degenerate_sizes() {
    assert!(EigenSolver::tridiag_eigenvalues(&[], &[]).is_empty());
    assert_eq!(EigenSolver::tridiag_eigenvalues(&[4.5], &[]), vec![4.5]);
    // Diagonal matrix: eigenvalues are the sorted diagonal.
    let ev = EigenSolver::tridiag_eigenvalues(&[3.0, -1.0, 2.0], &[0.0, 0.0]);
    assert_eq!(ev, vec![-1.0, 2.0, 3.0]);
}

#[test]
fn tridiag_eigenvalues_match_closed_form_for_larger_laplacians() {
    for n in [8_usize, 25, 60] {
        let alpha = vec![2.0; n];
        let beta = vec![-1.0; n - 1];
        let ev = EigenSolver::tridiag_eigenvalues(&alpha, &beta);
        assert_eq!(ev.len(), n);
        for (j, got) in ev.iter().enumerate() {
            let want = laplacian_eigenvalue(j + 1, n);
            assert!((got - want).abs() < 1e-12, "n = {n}, j = {}", j + 1);
        }
    }
}

// ---------------------------------------------------------------------------
// Hessenberg eigenvalues with imaginary parts
// ---------------------------------------------------------------------------

#[test]
fn hessenberg_complex_eigenvalues_of_a_rotation() {
    // [[1, -2], [2, 1]] has eigenvalues 1 +- 2i.
    let ev = EigenSolver::hessenberg_eigenvalues_complex(&[1.0, -2.0, 2.0, 1.0], 2);
    assert_eq!(ev.len(), 2);
    assert!((ev[0].0 - 1.0).abs() < 1e-12 && (ev[0].1 + 2.0).abs() < 1e-12);
    assert!((ev[1].0 - 1.0).abs() < 1e-12 && (ev[1].1 - 2.0).abs() < 1e-12);
}

#[test]
fn hessenberg_complex_eigenvalues_of_the_cube_roots_of_unity() {
    // Companion matrix of x^3 - 1: eigenvalues 1 and -1/2 +- i sqrt(3)/2.
    let h = [0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
    let ev = EigenSolver::hessenberg_eigenvalues_complex(&h, 3);
    let half_sqrt3 = 0.75_f64.sqrt();
    let expected = [(-0.5, -half_sqrt3), (-0.5, half_sqrt3), (1.0, 0.0)];
    assert_eq!(ev.len(), 3);
    for (got, want) in ev.iter().zip(expected.iter()) {
        assert!(
            (got.0 - want.0).abs() < 1e-10 && (got.1 - want.1).abs() < 1e-10,
            "got {got:?}, want {want:?}"
        );
    }
    // The real-part wrapper keeps working on the same input.
    let real_parts = EigenSolver::hessenberg_eigenvalues(&h, 3);
    assert!((real_parts[0] + 0.5).abs() < 1e-10);
    assert!((real_parts[1] + 0.5).abs() < 1e-10);
    assert!((real_parts[2] - 1.0).abs() < 1e-10);
}

#[test]
fn hessenberg_complex_eigenvalues_of_a_block_triangular_matrix() {
    // Diagonal blocks [[2, -3], [3, 2]] (eigenvalues 2 +- 3i), [1] and [5].
    let h = [
        2.0, -3.0, 1.0, 0.5, //
        3.0, 2.0, 0.0, 2.0, //
        0.0, 0.0, 1.0, 4.0, //
        0.0, 0.0, 0.0, 5.0,
    ];
    let ev = EigenSolver::hessenberg_eigenvalues_complex(&h, 4);
    let expected = [(1.0, 0.0), (2.0, -3.0), (2.0, 3.0), (5.0, 0.0)];
    assert_eq!(ev.len(), 4);
    for (got, want) in ev.iter().zip(expected.iter()) {
        assert!(
            (got.0 - want.0).abs() < 1e-12 && (got.1 - want.1).abs() < 1e-12,
            "got {got:?}, want {want:?}"
        );
    }
}

#[test]
fn hessenberg_eigenvalues_trivial_sizes() {
    assert!(EigenSolver::hessenberg_eigenvalues_complex(&[], 0).is_empty());
    assert_eq!(
        EigenSolver::hessenberg_eigenvalues_complex(&[-2.5], 1),
        vec![(-2.5, 0.0)]
    );
    // A zero matrix deflates immediately.
    let zeros = EigenSolver::hessenberg_eigenvalues_complex(&[0.0; 9], 3);
    assert!(zeros.iter().all(|e| *e == (0.0, 0.0)));
}

// ---------------------------------------------------------------------------
// F0790: Arnoldi / Lanczos
// ---------------------------------------------------------------------------

#[test]
fn lanczos_two_largest_of_laplacian_30_match_closed_form() {
    let n = 30;
    let lap = laplacian(n);
    let res = EigenSolver::lanczos(&lap, 2, n, 1e-10);
    assert!(res.converged, "Lanczos did not converge");
    assert_eq!(res.eigenvalues.len(), 2);
    assert_eq!(res.eigenvectors.len(), 2);
    // Ascending: the second largest first.  (j = 30 is an antisymmetric mode that
    // a constant start vector could never reach.)
    let expected = [laplacian_eigenvalue(29, n), laplacian_eigenvalue(30, n)];
    for (i, want) in expected.iter().enumerate() {
        let (got, v) = (res.eigenvalues[i], &res.eigenvectors[i]);
        assert!((got - want).abs() < 1e-9, "eigenvalue {i}: {got} vs {want}");
        assert!((norm(v) - 1.0).abs() < 1e-10, "Ritz vector {i} not unit");
        assert!(
            residual_norm(&lap, got, v) < 1e-6,
            "residual {i} = {}",
            residual_norm(&lap, got, v)
        );
    }
    assert!(inner(&res.eigenvectors[0], &res.eigenvectors[1]).abs() < 1e-8);
}

#[test]
fn arnoldi_two_largest_of_laplacian_30_match_closed_form() {
    let n = 30;
    let lap = laplacian(n);
    let res = EigenSolver::arnoldi(&lap, 2, n, 1e-10);
    assert!(res.converged, "Arnoldi did not converge");
    assert_eq!(res.eigenvalues.len(), 2);
    let expected = [laplacian_eigenvalue(29, n), laplacian_eigenvalue(30, n)];
    for (i, want) in expected.iter().enumerate() {
        let got = res.eigenvalues[i];
        assert!((got - want).abs() < 1e-9, "eigenvalue {i}: {got} vs {want}");
        assert!(residual_norm(&lap, got, &res.eigenvectors[i]) < 1e-6);
    }
}

#[test]
fn lanczos_full_krylov_space_recovers_every_eigenvalue() {
    let n = 8;
    let lap = laplacian(n);
    let res = EigenSolver::lanczos(&lap, 8, 8, 1e-12);
    assert!(res.converged);
    assert_eq!(res.eigenvalues.len(), 8);
    assert_eq!(res.iterations, 8);
    for (j, got) in res.eigenvalues.iter().enumerate() {
        let want = laplacian_eigenvalue(j + 1, n);
        assert!((got - want).abs() < 1e-8, "j = {}: {got} vs {want}", j + 1);
        assert!(residual_norm(&lap, *got, &res.eigenvectors[j]) < 1e-10);
    }
    // The eigenvectors of distinct eigenvalues are orthonormal.
    for i in 0..8 {
        for j in 0..8 {
            let expected = if i == j { 1.0 } else { 0.0 };
            let ip = inner(&res.eigenvectors[i], &res.eigenvectors[j]);
            assert!((ip - expected).abs() < 1e-9, "<v{i}, v{j}> = {ip}");
        }
    }
}

#[test]
fn arnoldi_full_krylov_space_recovers_every_eigenvalue() {
    let n = 8;
    let lap = laplacian(n);
    let res = EigenSolver::arnoldi(&lap, 8, 8, 1e-12);
    assert!(res.converged);
    assert_eq!(res.eigenvalues.len(), 8);
    for (j, got) in res.eigenvalues.iter().enumerate() {
        let want = laplacian_eigenvalue(j + 1, n);
        assert!((got - want).abs() < 1e-8, "j = {}: {got} vs {want}", j + 1);
        assert!(residual_norm(&lap, *got, &res.eigenvectors[j]) < 1e-8);
    }
}

#[test]
fn lanczos_start_vector_reaches_antisymmetric_modes() {
    // The largest eigenvalue of the 10-point Laplacian belongs to the
    // antisymmetric mode j = 10; a constant start vector has no component along
    // it and would report the j = 9 eigenvalue instead.
    let lap = laplacian(10);
    let res = EigenSolver::lanczos(&lap, 1, 10, 1e-10);
    assert!(res.converged);
    assert!((res.eigenvalues[0] - laplacian_eigenvalue(10, 10)).abs() < 1e-9);
}

#[test]
fn lanczos_and_arnoldi_report_non_convergence_honestly() {
    let n = 100;
    let lap = laplacian(n);
    for res in [
        EigenSolver::lanczos(&lap, 3, 7, 1e-10),
        EigenSolver::arnoldi(&lap, 3, 7, 1e-10),
    ] {
        assert!(
            !res.converged,
            "7 Krylov steps cannot resolve three eigenpairs"
        );
        assert_eq!(res.eigenvalues.len(), 3);
        assert_eq!(res.iterations, 7);
        // Ritz values stay inside the spectrum and the residuals are genuinely large.
        let lo = laplacian_eigenvalue(1, n);
        let hi = laplacian_eigenvalue(n, n);
        for (lambda, v) in res.eigenvalues.iter().zip(res.eigenvectors.iter()) {
            assert!(*lambda > lo && *lambda < hi + 1e-12);
            assert!((norm(v) - 1.0).abs() < 1e-10);
        }
        let worst = res
            .eigenvalues
            .iter()
            .zip(res.eigenvectors.iter())
            .map(|(l, v)| residual_norm(&lap, *l, v))
            .fold(0.0_f64, f64::max);
        assert!(worst > 1e-6, "worst residual {worst}");
    }
}

#[test]
fn lanczos_finds_repeated_eigenvalues_by_restarting() {
    // A single-vector Krylov space sees each distinct eigenvalue once; the
    // restart after the invariant subspace is exhausted must find the copies.
    let entries = [1.0, 1.0, 1.0, 2.0, 2.0, 5.0, 5.0, 5.0, 5.0, 7.0];
    let a = diagonal_matrix(&entries);
    let res = EigenSolver::lanczos(&a, 4, 10, 1e-10);
    assert!(res.converged);
    let expected = [5.0, 5.0, 5.0, 7.0];
    assert_eq!(res.eigenvalues.len(), 4);
    for (got, want) in res.eigenvalues.iter().zip(expected.iter()) {
        assert!((got - want).abs() < 1e-10, "got {got}, want {want}");
    }
    for i in 0..4 {
        assert!(residual_norm(&a, res.eigenvalues[i], &res.eigenvectors[i]) < 1e-10);
        for j in 0..i {
            let ip = inner(&res.eigenvectors[i], &res.eigenvectors[j]);
            assert!(ip.abs() < 1e-9, "<v{i}, v{j}> = {ip}");
        }
    }
}

#[test]
fn arnoldi_finds_repeated_eigenvalues_with_independent_vectors() {
    let entries = [1.0, 1.0, 1.0, 2.0, 2.0, 5.0, 5.0, 5.0, 5.0, 7.0];
    let a = diagonal_matrix(&entries);
    let res = EigenSolver::arnoldi(&a, 4, 10, 1e-10);
    assert!(res.converged);
    let expected = [5.0, 5.0, 5.0, 7.0];
    for (got, want) in res.eigenvalues.iter().zip(expected.iter()) {
        assert!((got - want).abs() < 1e-10, "got {got}, want {want}");
    }
    // The three copies of the eigenvalue 5 must not collapse onto one vector.
    for i in 0..3 {
        for j in 0..i {
            let ip = inner(&res.eigenvectors[i], &res.eigenvectors[j]);
            assert!(ip.abs() < 1e-8, "<v{i}, v{j}> = {ip}");
        }
    }
    for i in 0..4 {
        assert!(residual_norm(&a, res.eigenvalues[i], &res.eigenvectors[i]) < 1e-10);
    }
}

#[test]
fn arnoldi_non_symmetric_dominant_eigenpairs() {
    let n = 20;
    let c = 0.3;
    let a = convection_diffusion(n, c);
    let res = EigenSolver::arnoldi(&a, 3, n, 1e-9);
    assert!(res.converged);
    let expected = [
        convection_diffusion_eigenvalue(3, n, c),
        convection_diffusion_eigenvalue(2, n, c),
        convection_diffusion_eigenvalue(1, n, c),
    ];
    assert_eq!(res.eigenvalues.len(), 3);
    for (i, want) in expected.iter().enumerate() {
        let got = res.eigenvalues[i];
        assert!((got - want).abs() < 1e-8, "eigenvalue {i}: {got} vs {want}");
        assert!(residual_norm(&a, got, &res.eigenvectors[i]) < 1e-8);
    }
}

#[test]
fn arnoldi_skips_complex_conjugate_pairs() {
    // Eigenvalues: +-i (a rotation block), 3 and 0.5.  Complex pairs have no real
    // Ritz vector, so only the real eigenpairs can be reported.
    let a = from_dense(&[
        vec![0.0, -1.0, 0.0, 0.0],
        vec![1.0, 0.0, 0.0, 0.0],
        vec![0.0, 0.0, 3.0, 0.0],
        vec![0.0, 0.0, 0.0, 0.5],
    ]);
    let res = EigenSolver::arnoldi(&a, 2, 4, 1e-10);
    assert!(res.converged);
    assert_eq!(res.eigenvalues.len(), 2);
    assert!((res.eigenvalues[0] - 0.5).abs() < 1e-10);
    assert!((res.eigenvalues[1] - 3.0).abs() < 1e-10);
    for i in 0..2 {
        assert!(residual_norm(&a, res.eigenvalues[i], &res.eigenvectors[i]) < 1e-10);
    }
    // Asking for three real pairs when only two exist cannot be satisfied.
    let res3 = EigenSolver::arnoldi(&a, 3, 4, 1e-10);
    assert!(!res3.converged);
    assert_eq!(res3.eigenvalues.len(), 2);
}

#[test]
fn eigensolvers_handle_trivial_requests() {
    let lap = laplacian(6);
    for res in [
        EigenSolver::lanczos(&lap, 0, 10, 1e-8),
        EigenSolver::arnoldi(&lap, 0, 10, 1e-8),
        EigenSolver::lobpcg(&lap, 0, 1e-8, 10),
    ] {
        assert!(res.eigenvalues.is_empty() && res.eigenvectors.is_empty());
        assert!(res.converged);
    }
    let empty = SparseCSR::new(0, 0);
    assert!(
        EigenSolver::lanczos(&empty, 3, 10, 1e-8)
            .eigenvalues
            .is_empty()
    );
    assert!(
        EigenSolver::arnoldi(&empty, 3, 10, 1e-8)
            .eigenvalues
            .is_empty()
    );
    assert!(
        EigenSolver::lobpcg(&empty, 3, 1e-8, 10)
            .eigenvalues
            .is_empty()
    );
    // k larger than n is clamped.
    let res = EigenSolver::lanczos(&laplacian(3), 10, 10, 1e-10);
    assert_eq!(res.eigenvalues.len(), 3);
    assert!(res.converged);
    // 1 x 1 matrix.
    let single = diagonal_matrix(&[4.0]);
    let res = EigenSolver::lanczos(&single, 1, 5, 1e-12);
    assert!((res.eigenvalues[0] - 4.0).abs() < 1e-12);
    let res = EigenSolver::arnoldi(&single, 1, 5, 1e-12);
    assert!((res.eigenvalues[0] - 4.0).abs() < 1e-12);
}

#[test]
fn eigensolvers_reject_non_square_matrices_without_panicking() {
    let rectangular = SparseCSR::from_triplets(3, 5, &[0, 1, 2], &[0, 3, 4], &[1.0, 2.0, 3.0]);
    for res in [
        EigenSolver::lanczos(&rectangular, 2, 10, 1e-8),
        EigenSolver::arnoldi(&rectangular, 2, 10, 1e-8),
        EigenSolver::lobpcg(&rectangular, 2, 1e-8, 10),
    ] {
        assert!(res.eigenvalues.is_empty() && res.eigenvectors.is_empty());
        assert!(!res.converged);
    }
}

#[test]
fn lanczos_zero_matrix_yields_zero_eigenvalues() {
    let a = SparseCSR::new(5, 5);
    let res = EigenSolver::lanczos(&a, 2, 5, 1e-10);
    assert!(res.converged);
    assert_eq!(res.eigenvalues.len(), 2);
    assert!(res.eigenvalues.iter().all(|v| v.abs() < 1e-12));
    for i in 0..2 {
        assert!((norm(&res.eigenvectors[i]) - 1.0).abs() < 1e-10);
    }
}

// ---------------------------------------------------------------------------
// F0792: LOBPCG
// ---------------------------------------------------------------------------

#[test]
fn lobpcg_scaled_laplacian_three_smallest() {
    // 10 * laplacian_1d(50): the old fixed-step gradient descent reported
    // converged = true with values from the middle of the spectrum once the norm
    // exceeded 20.
    let n = 50;
    let a = scaled(&laplacian(n), 10.0);
    let res = EigenSolver::lobpcg(&a, 3, 1e-8, 1000);
    assert!(res.converged, "LOBPCG did not converge");
    assert_eq!(res.eigenvalues.len(), 3);
    assert_eq!(res.eigenvectors.len(), 3);
    for j in 1..=3 {
        let want = 10.0 * laplacian_eigenvalue(j, n);
        let got = res.eigenvalues[j - 1];
        assert!((got - want).abs() < 1e-6, "j = {j}: {got} vs {want}");
        let v = &res.eigenvectors[j - 1];
        assert!((norm(v) - 1.0).abs() < 1e-10);
        assert!(residual_norm(&a, got, v) < 1e-6);
    }
}

#[test]
fn lobpcg_unscaled_laplacian_two_smallest_are_accurate() {
    let n = 80;
    let a = laplacian(n);
    let res = EigenSolver::lobpcg(&a, 2, 1e-9, 2000);
    assert!(res.converged);
    for j in 1..=2 {
        let want = laplacian_eigenvalue(j, n);
        assert!((res.eigenvalues[j - 1] - want).abs() < 1e-8, "j = {j}");
    }
}

#[test]
fn lobpcg_is_not_fooled_by_a_diagonal_matrix() {
    // Unit vectors are exact eigenvectors of a diagonal matrix, so a unit-vector
    // initial block would be a fixed point with the wrong eigenvalues.
    let entries: Vec<f64> = (0..60).map(|i| 100.0 - i as f64).collect();
    let a = diagonal_matrix(&entries);
    let res = EigenSolver::lobpcg(&a, 3, 1e-9, 500);
    assert!(res.converged);
    let expected = [41.0, 42.0, 43.0];
    for (got, want) in res.eigenvalues.iter().zip(expected.iter()) {
        assert!((got - want).abs() < 1e-8, "got {got}, want {want}");
    }
}

#[test]
fn lobpcg_two_dimensional_grid_with_degenerate_eigenvalues() {
    // 8 x 8 grid Laplacian: eigenvalues are sums of 1-D eigenvalues and come in
    // degenerate pairs.  Reference: dense symmetric solver on the same matrix.
    let m = 8;
    let n = m * m;
    let (mut r, mut c, mut v) = (Vec::new(), Vec::new(), Vec::new());
    let mut dense = vec![vec![0.0_f64; n]; n];
    for i in 0..m {
        for j in 0..m {
            let p = i * m + j;
            r.push(p);
            c.push(p);
            v.push(4.0);
            dense[p][p] = 4.0;
            let neighbours = [
                (i > 0).then(|| p - m),
                (i + 1 < m).then(|| p + m),
                (j > 0).then(|| p - 1),
                (j + 1 < m).then(|| p + 1),
            ];
            for q in neighbours.into_iter().flatten() {
                r.push(p);
                c.push(q);
                v.push(-1.0);
                dense[p][q] = -1.0;
            }
        }
    }
    let a = SparseCSR::from_triplets(n, n, &r, &c, &v);
    let (reference, _) = symmetric_eigen_n(&dense);
    let res = EigenSolver::lobpcg(&a, 4, 1e-9, 1000);
    assert!(res.converged);
    for (i, want) in reference.iter().take(4).enumerate() {
        let got = res.eigenvalues[i];
        assert!((got - want).abs() < 1e-7, "eigenvalue {i}: {got} vs {want}");
        assert!(residual_norm(&a, got, &res.eigenvectors[i]) < 1e-6);
    }
}

#[test]
fn lobpcg_handles_an_indefinite_matrix() {
    // Laplacian - 2 I has eigenvalues -2 cos(j pi / (n + 1)): the smallest are negative.
    let n = 70;
    let shifted = {
        let lap = laplacian(n);
        let mut rows = Vec::new();
        let mut cols = Vec::new();
        let mut vals = Vec::new();
        for (i, bounds) in lap.row_ptr.windows(2).enumerate() {
            for k in bounds[0]..bounds[1] {
                rows.push(i);
                cols.push(lap.col_idx[k]);
                vals.push(lap.values[k]);
            }
            rows.push(i);
            cols.push(i);
            vals.push(-2.0);
        }
        SparseCSR::from_triplets(n, n, &rows, &cols, &vals)
    };
    let res = EigenSolver::lobpcg(&shifted, 2, 1e-9, 2000);
    assert!(res.converged);
    for j in 1..=2 {
        let want = laplacian_eigenvalue(j, n) - 2.0;
        assert!((res.eigenvalues[j - 1] - want).abs() < 1e-8, "j = {j}");
    }
}

#[test]
fn lobpcg_reports_non_convergence_and_ritz_values_bound_the_spectrum() {
    let n = 120;
    let a = laplacian(n);
    let res = EigenSolver::lobpcg(&a, 2, 1e-12, 3);
    assert!(!res.converged);
    assert_eq!(res.iterations, 3);
    assert_eq!(res.eigenvalues.len(), 2);
    // Rayleigh-Ritz values are upper bounds of the corresponding eigenvalues.
    for j in 1..=2 {
        assert!(res.eigenvalues[j - 1] >= laplacian_eigenvalue(j, n) - 1e-12);
    }
    assert!(res.eigenvalues[0] <= res.eigenvalues[1]);
}

#[test]
fn lobpcg_small_problem_uses_the_dense_solver() {
    let n = 6;
    let a = laplacian(n);
    let res = EigenSolver::lobpcg(&a, 2, 1e-12, 100);
    assert!(res.converged);
    assert_eq!(res.eigenvalues.len(), 2);
    for j in 1..=2 {
        assert!((res.eigenvalues[j - 1] - laplacian_eigenvalue(j, n)).abs() < 1e-12);
        assert!(residual_norm(&a, res.eigenvalues[j - 1], &res.eigenvectors[j - 1]) < 1e-12);
    }
}

// ---------------------------------------------------------------------------
// F0798: TFQMR
// ---------------------------------------------------------------------------

#[test]
fn tfqmr_laplacian_5_converges_with_a_true_residual() {
    // The old implementation returned converged = true with residual 1.09e19.
    let a = laplacian(5);
    let b: Vec<f64> = (1..=5).map(f64::from).collect();
    let res = IterativeSolver::tfqmr(&a, &b, None, 1e-10, 100);
    assert!(res.converged);
    assert!(res.residual < 1e-10, "residual {}", res.residual);
    assert!((true_residual(&a, &res.x, &b) - res.residual).abs() < 1e-12);
    // Agrees with the CG solution.
    let reference = IterativeSolver::conjugate_gradient(&a, &b, None, 1e-12, 100);
    for (u, w) in res.x.iter().zip(reference.x.iter()) {
        assert!((u - w).abs() < 1e-8, "{u} vs {w}");
    }
}

#[test]
fn tfqmr_laplacian_10_converges() {
    let a = laplacian(10);
    let b: Vec<f64> = (1..=10).map(f64::from).collect();
    let res = IterativeSolver::tfqmr(&a, &b, None, 1e-9, 200);
    assert!(res.converged);
    assert!(true_residual(&a, &res.x, &b) < 1e-9);
}

#[test]
fn tfqmr_solves_a_non_symmetric_convection_diffusion_system() {
    // Moderately non-normal (eigenvector condition number about 3e3): TFQMR
    // reaches the tolerance within n steps.  A strongly non-normal variant (c =
    // 0.4 gives about 2e5) makes the quasi-residual recurrences stall near 1e-5
    // in an independent double-precision port of Saad's Algorithm 7.8 as well,
    // which is why the next test only checks that the flag tells the truth.
    let n = 40;
    let a = convection_diffusion(n, 0.2);
    let x_true: Vec<f64> = (0..n).map(|i| 1.0 + (i as f64 * 0.3).sin()).collect();
    let b = a.matvec(&x_true);
    let res = IterativeSolver::tfqmr(&a, &b, None, 1e-10, 500);
    assert!(res.converged, "residual {}", res.residual);
    assert!(true_residual(&a, &res.x, &b) < 1e-10);
    for (u, w) in res.x.iter().zip(x_true.iter()) {
        assert!((u - w).abs() < 1e-6, "{u} vs {w}");
    }
}

#[test]
fn tfqmr_converged_flag_always_reflects_the_true_residual() {
    let a = convection_diffusion(30, 0.4);
    let b = vec![1.0; 30];
    for max_iter in [1_usize, 2, 5, 400] {
        let tol = 1e-8;
        let res = IterativeSolver::tfqmr(&a, &b, None, tol, max_iter);
        let actual = true_residual(&a, &res.x, &b);
        assert!((actual - res.residual).abs() <= 1e-10 * (1.0 + actual));
        assert_eq!(res.converged, actual < tol, "max_iter = {max_iter}");
    }
}

#[test]
fn tfqmr_trivial_right_hand_sides_and_initial_guesses() {
    let a = laplacian(5);
    let zero = vec![0.0; 5];
    let res = IterativeSolver::tfqmr(&a, &zero, None, 1e-10, 50);
    assert!(res.converged);
    assert_eq!(res.iterations, 0);
    assert!(res.x.iter().all(|v| *v == 0.0));
    let x_exact = vec![1.0, -2.0, 3.0, -4.0, 5.0];
    let b = a.matvec(&x_exact);
    let res = IterativeSolver::tfqmr(&a, &b, Some(&x_exact), 1e-10, 50);
    assert!(res.converged);
    assert_eq!(res.iterations, 0);
    // A decent initial guess is honoured and improved on.
    let guess: Vec<f64> = x_exact.iter().map(|v| v + 0.01).collect();
    let res = IterativeSolver::tfqmr(&a, &b, Some(&guess), 1e-10, 50);
    assert!(res.converged);
    assert!(true_residual(&a, &res.x, &b) < 1e-10);
}

#[test]
fn tfqmr_breakdown_is_reported_not_hidden() {
    // A = 0 makes r~^T v vanish immediately; the solver must stop and say so.
    let a = SparseCSR::new(3, 3);
    let b = vec![1.0, 2.0, 2.0];
    let res = IterativeSolver::tfqmr(&a, &b, None, 1e-8, 20);
    assert!(!res.converged);
    assert!((res.residual - 3.0).abs() < 1e-12);
    assert!(res.x.iter().all(|v| v.is_finite()));
}

// ---------------------------------------------------------------------------
// F0802: matrix logarithm
// ---------------------------------------------------------------------------

#[test]
fn logm_of_a_diagonal_matrix_is_the_diagonal_of_logarithms() {
    // The old Padé coefficients returned about half the logarithm.
    let log_a = MatrixFunctions::logm(&[4.0, 0.0, 0.0, 9.0], 2);
    assert!((log_a[0] - 4.0_f64.ln()).abs() < 1e-9, "{}", log_a[0]);
    assert!((log_a[3] - 9.0_f64.ln()).abs() < 1e-9, "{}", log_a[3]);
    assert!(log_a[1].abs() < 1e-9 && log_a[2].abs() < 1e-9);
}

#[test]
fn logm_of_scalars() {
    assert!((MatrixFunctions::logm(&[1.1], 1)[0] - 1.1_f64.ln()).abs() < 1e-15);
    assert!((MatrixFunctions::logm(&[4.0], 1)[0] - 4.0_f64.ln()).abs() < 1e-10);
    assert!((MatrixFunctions::logm(&[0.2], 1)[0] - 0.2_f64.ln()).abs() < 1e-10);
    assert!(MatrixFunctions::logm(&[1.0], 1)[0].abs() < 1e-15);
}

#[test]
fn logm_of_identity_and_empty() {
    let log_i = MatrixFunctions::logm(&MatrixFunctions::eye(3), 3);
    assert!(log_i.iter().all(|v| v.abs() < 1e-14));
    assert!(MatrixFunctions::logm(&[], 0).is_empty());
}

#[test]
fn logm_is_the_inverse_of_expm_for_a_spd_matrix() {
    let a = [4.0, 1.0, 0.5, 1.0, 3.0, 0.2, 0.5, 0.2, 2.0];
    let log_a = MatrixFunctions::logm(&a, 3);
    let back = MatrixFunctions::expm(&log_a, 3);
    for (got, want) in back.iter().zip(a.iter()) {
        assert!((got - want).abs() < 1e-8, "{got} vs {want}");
    }
    // The logarithm of a symmetric matrix is symmetric.
    for i in 0..3 {
        for j in 0..3 {
            assert!((log_a[i * 3 + j] - log_a[j * 3 + i]).abs() < 1e-10);
        }
    }
}

#[test]
fn logm_inverts_expm_for_a_symmetric_argument() {
    let b = [0.3, -0.2, 0.1, -0.2, 0.5, 0.25, 0.1, 0.25, -0.4];
    let exp_b = MatrixFunctions::expm(&b, 3);
    let log_exp_b = MatrixFunctions::logm(&exp_b, 3);
    for (got, want) in log_exp_b.iter().zip(b.iter()) {
        assert!((got - want).abs() < 1e-9, "{got} vs {want}");
    }
}

#[test]
fn logm_of_a_rotated_diagonal_matches_the_closed_form() {
    // A = Q diag(2, 5) Q^T  =>  log A = Q diag(ln 2, ln 5) Q^T.
    let angle = 0.7_f64;
    let (c, s) = (angle.cos(), angle.sin());
    let (l1, l2) = (2.0_f64, 5.0_f64);
    let build = |d1: f64, d2: f64| {
        [
            c * c * d1 + s * s * d2,
            c * s * (d1 - d2),
            c * s * (d1 - d2),
            s * s * d1 + c * c * d2,
        ]
    };
    let a = build(l1, l2);
    let expected = build(l1.ln(), l2.ln());
    let log_a = MatrixFunctions::logm(&a, 2);
    for (got, want) in log_a.iter().zip(expected.iter()) {
        assert!((got - want).abs() < 1e-10, "{got} vs {want}");
    }
}

#[test]
fn logm_handles_widely_scaled_spectra() {
    let log_a = MatrixFunctions::logm(&[1.0e4, 0.0, 0.0, 1.0e-4], 2);
    assert!((log_a[0] - 1.0e4_f64.ln()).abs() < 1e-8, "{}", log_a[0]);
    assert!((log_a[3] - 1.0e-4_f64.ln()).abs() < 1e-8, "{}", log_a[3]);
}

#[test]
fn logm_of_a_scalar_multiple_shifts_by_the_log_of_the_scalar() {
    let a = [3.0, 0.5, 0.5, 2.0];
    let log_a = MatrixFunctions::logm(&a, 2);
    let scaled_a: Vec<f64> = a.iter().map(|v| 6.0 * v).collect();
    let log_scaled = MatrixFunctions::logm(&scaled_a, 2);
    let shift = 6.0_f64.ln();
    for i in 0..2 {
        for j in 0..2 {
            let expected = log_a[i * 2 + j] + if i == j { shift } else { 0.0 };
            assert!((log_scaled[i * 2 + j] - expected).abs() < 1e-9);
        }
    }
}

// ---------------------------------------------------------------------------
// expm: the [13/13] Pade coefficients used by logm round trips
// ---------------------------------------------------------------------------

#[test]
fn expm_of_the_identity_and_diagonals_is_accurate_to_rounding() {
    // A hand-typed coefficient table used to limit exp(I) to about 7 digits.
    let e = MatrixFunctions::expm(&MatrixFunctions::eye(3), 3);
    for i in 0..3 {
        assert!(
            (e[i * 3 + i] - std::f64::consts::E).abs() < 1e-14,
            "exp(1) = {}",
            e[i * 3 + i]
        );
    }
    let d = [0.5, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 3.0];
    let ed = MatrixFunctions::expm(&d, 3);
    let expected = [0.5_f64.exp(), (-1.0_f64).exp(), 3.0_f64.exp()];
    for i in 0..3 {
        assert!(
            (ed[i * 3 + i] - expected[i]).abs() < 1e-13 * expected[i].max(1.0),
            "exp({}) = {}",
            d[i * 3 + i],
            ed[i * 3 + i]
        );
    }
}

#[test]
fn expm_of_a_rotation_generator_is_a_rotation() {
    let theta = 1.3_f64;
    let generator = [0.0, -theta, theta, 0.0];
    let r = MatrixFunctions::expm(&generator, 2);
    let expected = [theta.cos(), -theta.sin(), theta.sin(), theta.cos()];
    for (got, want) in r.iter().zip(expected.iter()) {
        assert!((got - want).abs() < 1e-14, "{got} vs {want}");
    }
}

// ---------------------------------------------------------------------------
// F0742: MINRES
// ---------------------------------------------------------------------------

/// Dense symmetric matrix (all rows given) as a CSR matrix.
fn symmetric_from_dense(rows: &[Vec<f64>]) -> SparseCSR {
    from_dense(rows)
}

#[test]
fn minres_solves_a_one_by_one_system_exactly() {
    // A = [2], b = [2]: the old recurrence returned b / (1 + b^2) = 0.4.
    let a = diagonal_matrix(&[2.0]);
    let res = IterativeSolver::minres(&a, &[2.0], None, 1e-12, 10);
    assert!(res.converged);
    assert!((res.x[0] - 1.0).abs() < 1e-14, "x = {}", res.x[0]);
    assert!(res.residual < 1e-14);
    assert_eq!(res.iterations, 1);
}

#[test]
fn minres_matches_conjugate_gradient_on_the_laplacian() {
    let n = 20;
    let a = laplacian(n);
    let b: Vec<f64> = (0..n).map(|i| (i + 1) as f64).collect();
    let cg = IterativeSolver::conjugate_gradient(&a, &b, None, 1e-12, 500);
    let mr = IterativeSolver::minres(&a, &b, None, 1e-9, 500);
    assert!(cg.converged && mr.converged);
    for i in 0..n {
        assert!(
            (cg.x[i] - mr.x[i]).abs() < 1e-8,
            "component {i}: cg {} vs minres {}",
            cg.x[i],
            mr.x[i]
        );
    }
    assert!(true_residual(&a, &mr.x, &b) < 1e-9);
}

#[test]
fn minres_solves_an_indefinite_diagonal_system() {
    // diag(-1, 2, 3) x = (1, 2, 3)  =>  x = (-1, 1, 1).
    let a = diagonal_matrix(&[-1.0, 2.0, 3.0]);
    let res = IterativeSolver::minres(&a, &[1.0, 2.0, 3.0], None, 1e-12, 50);
    assert!(res.converged, "residual {}", res.residual);
    let expected = [-1.0, 1.0, 1.0];
    for (got, want) in res.x.iter().zip(expected.iter()) {
        assert!((got - want).abs() < 1e-12, "got {got}, want {want}");
    }
    assert!(res.iterations <= 3, "three distinct eigenvalues");
}

#[test]
fn minres_solves_a_shifted_indefinite_laplacian() {
    // L - 2.3 I has eigenvalues of both signs; CG is not applicable.
    let n = 25;
    let mut rows = vec![vec![0.0; n]; n];
    for i in 0..n {
        rows[i][i] = 2.0 - 2.3;
        if i > 0 {
            rows[i][i - 1] = -1.0;
        }
        if i + 1 < n {
            rows[i][i + 1] = -1.0;
        }
    }
    let a = symmetric_from_dense(&rows);
    let b: Vec<f64> = (0..n).map(|i| (i as f64).sin() + 0.5).collect();
    let res = IterativeSolver::minres(&a, &b, None, 1e-9, 500);
    assert!(res.converged, "residual {}", res.residual);
    // Independent reference: dense Gaussian elimination.
    let reference = crate::numerical_methods::lu_solve(&rows, &b);
    for (i, (got, want)) in res.x.iter().zip(reference.iter()).enumerate() {
        assert!(
            (got - want).abs() < 1e-7 * (1.0 + want.abs()),
            "component {i}: {got} vs {want}"
        );
    }
}

#[test]
fn minres_reports_the_true_residual() {
    let a = laplacian(30);
    let b: Vec<f64> = (0..30).map(|i| ((i * 7 % 11) as f64) - 5.0).collect();
    for max_iter in [3usize, 12, 60] {
        let res = IterativeSolver::minres(&a, &b, None, 1e-10, max_iter);
        let actual = true_residual(&a, &res.x, &b);
        assert!(
            (res.residual - actual).abs() <= 1e-12 * (1.0 + actual),
            "max_iter {max_iter}: reported {} vs actual {actual}",
            res.residual
        );
        assert_eq!(res.converged, actual < 1e-10, "max_iter {max_iter}");
    }
}

#[test]
fn minres_residual_never_increases_with_the_iteration_count() {
    let a = laplacian(30);
    let b: Vec<f64> = (0..30).map(|i| ((i * 7 % 11) as f64) - 5.0).collect();
    let mut previous = f64::INFINITY;
    for k in 1..=30usize {
        let res = IterativeSolver::minres(&a, &b, None, 1e-300, k);
        assert_eq!(res.iterations, k);
        assert!(
            res.residual <= previous + 1e-12,
            "residual grew at iteration {k}: {previous} -> {}",
            res.residual
        );
        previous = res.residual;
    }
    assert!(previous < 1e-8, "30 steps solve a 30x30 system: {previous}");
}

#[test]
fn minres_honours_a_nonzero_initial_guess() {
    let a = laplacian(10);
    let b = vec![1.0; 10];
    let x0 = vec![0.3; 10];
    let from_guess = IterativeSolver::minres(&a, &b, Some(&x0), 1e-10, 100);
    let from_zero = IterativeSolver::minres(&a, &b, None, 1e-10, 100);
    assert!(from_guess.converged && from_zero.converged);
    for i in 0..10 {
        assert!((from_guess.x[i] - from_zero.x[i]).abs() < 1e-8);
    }
}

#[test]
fn minres_returns_immediately_for_an_exact_start_and_for_a_zero_rhs() {
    let a = laplacian(8);
    let b: Vec<f64> = (0..8).map(|i| (i + 1) as f64).collect();
    let exact = IterativeSolver::minres(&a, &b, None, 1e-12, 100);
    let restarted = IterativeSolver::minres(&a, &b, Some(&exact.x), 1e-6, 100);
    assert_eq!(restarted.iterations, 0);
    assert!(restarted.converged);
    // A zero right-hand side is solved by x = 0 even with a zero tolerance.
    let zero = IterativeSolver::minres(&a, &[0.0; 8], None, 0.0, 100);
    assert!(zero.converged);
    assert_eq!(zero.iterations, 0);
    assert!(zero.x.iter().all(|&v| v == 0.0));
}

#[test]
fn minres_solves_a_singular_consistent_system() {
    // diag(1, 2, 0) x = (1, 2, 0) is consistent; the minimum-residual
    // solution reached from x0 = 0 is (1, 1, 0).
    let a = SparseCSR::from_triplets(3, 3, &[0, 1], &[0, 1], &[1.0, 2.0]);
    let res = IterativeSolver::minres(&a, &[1.0, 2.0, 0.0], None, 1e-12, 50);
    assert!(res.converged, "residual {}", res.residual);
    let expected = [1.0, 1.0, 0.0];
    for (got, want) in res.x.iter().zip(expected.iter()) {
        assert!((got - want).abs() < 1e-12, "got {got}, want {want}");
    }
}

#[test]
fn minres_finishes_in_two_steps_for_two_distinct_eigenvalues() {
    // Krylov dimension equals the number of distinct eigenvalues (here 2).
    let a = diagonal_matrix(&[2.0, 2.0, 2.0, 5.0, 5.0]);
    let res = IterativeSolver::minres(&a, &[1.0; 5], None, 1e-12, 50);
    assert!(res.converged);
    assert_eq!(res.iterations, 2);
    let expected = [0.5, 0.5, 0.5, 0.2, 0.2];
    for (got, want) in res.x.iter().zip(expected.iter()) {
        assert!((got - want).abs() < 1e-13, "got {got}, want {want}");
    }
}

#[test]
fn minres_reports_non_convergence_honestly() {
    let a = laplacian(50);
    let b = vec![1.0; 50];
    let res = IterativeSolver::minres(&a, &b, None, 1e-12, 3);
    assert!(!res.converged);
    assert_eq!(res.iterations, 3);
    assert!(res.residual > 1e-6);
    assert!((res.residual - true_residual(&a, &res.x, &b)).abs() < 1e-12);
}

#[test]
fn minres_is_not_fooled_by_a_huge_right_hand_side() {
    // ||r0|| is far larger than ||A||; the Lanczos breakdown test must be
    // relative to ||A||, not to ||r0||.
    let a = diagonal_matrix(&[1e-3, 2e-3]);
    let res = IterativeSolver::minres(&a, &[1e10, 1e10], None, 1e-3, 50);
    assert!(res.converged, "residual {}", res.residual);
    assert!((res.x[0] - 1e13).abs() < 1e13 * 1e-9);
    assert!((res.x[1] - 5e12).abs() < 5e12 * 1e-9);
}

// ---------------------------------------------------------------------------
// F0733: svd_small and its consumers
// ---------------------------------------------------------------------------

use super::{LowRankApproximation, TensorDecompResult, TensorDecomposition, qr_thin, svd_small};

/// Deterministic pseudo-random `m x n` matrix with entries in `[-1, 1]`.
fn pseudo_random_matrix(m: usize, n: usize, seed: u64) -> Vec<f64> {
    let mut state = seed;
    (0..m * n)
        .map(|_| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((state >> 11) as f64 / (1u64 << 53) as f64) * 2.0 - 1.0
        })
        .collect()
}

/// `max |(M^T M) - I|` for a row-major `rows x cols` matrix (orthonormal columns).
fn column_orthonormality_error(mat: &[f64], rows: usize, cols: usize) -> f64 {
    let mut worst = 0.0_f64;
    for p in 0..cols {
        for q in 0..cols {
            let d: f64 = (0..rows)
                .map(|i| mat[i * cols + p] * mat[i * cols + q])
                .sum();
            let target = if p == q { 1.0 } else { 0.0 };
            worst = worst.max((d - target).abs());
        }
    }
    worst
}

/// `max |(M M^T) - I|` for a row-major `rows x cols` matrix (orthonormal rows).
fn row_orthonormality_error(mat: &[f64], rows: usize, cols: usize) -> f64 {
    let mut worst = 0.0_f64;
    for p in 0..rows {
        for q in 0..rows {
            let d: f64 = (0..cols)
                .map(|j| mat[p * cols + j] * mat[q * cols + j])
                .sum();
            let target = if p == q { 1.0 } else { 0.0 };
            worst = worst.max((d - target).abs());
        }
    }
    worst
}

/// `max |U diag(s) Vt - A|` for rank-`k` factors of an `m x n` matrix
/// (`U` is `m x k`, `Vt` is `k x n`).
fn svd_reconstruction_error(
    a: &[f64],
    (m, n, k): (usize, usize, usize),
    u: &[f64],
    s: &[f64],
    vt: &[f64],
) -> f64 {
    let mut worst = 0.0_f64;
    for i in 0..m {
        for j in 0..n {
            let value: f64 = (0..k).map(|l| u[i * k + l] * s[l] * vt[l * n + j]).sum();
            worst = worst.max((value - a[i * n + j]).abs());
        }
    }
    worst
}

fn assert_valid_thin_svd(a: &[f64], m: usize, n: usize, tolerance: f64) {
    let k = m.min(n);
    let (u, s, vt) = svd_small(a, m, n);
    assert_eq!(u.len(), m * k, "U must be m x k");
    assert_eq!(s.len(), k);
    assert_eq!(vt.len(), k * n, "Vt must be k x n");
    assert!(
        s.iter().all(|&v| v >= 0.0),
        "singular values are nonnegative"
    );
    assert!(
        s.windows(2).all(|w| w[0] >= w[1]),
        "singular values must be sorted descending: {s:?}"
    );
    let reconstruction = svd_reconstruction_error(a, (m, n, k), &u, &s, &vt);
    assert!(
        reconstruction < tolerance,
        "{m}x{n}: reconstruction error {reconstruction}"
    );
    assert!(
        column_orthonormality_error(&u, m, k) < tolerance,
        "{m}x{n}: U columns are not orthonormal"
    );
    assert!(
        row_orthonormality_error(&vt, k, n) < tolerance,
        "{m}x{n}: Vt rows are not orthonormal"
    );
}

#[test]
fn svd_small_reconstructs_random_wide_and_tall_matrices() {
    // The old routine never rotated the columns beyond min(m, n) of a wide
    // matrix, left U/V unsorted and returned V where Vt was documented.
    for &(m, n, seed) in &[
        (3usize, 5usize, 1u64),
        (5, 3, 2),
        (4, 4, 3),
        (1, 6, 4),
        (6, 1, 5),
        (2, 9, 6),
        (9, 2, 7),
        (7, 8, 8),
        (12, 5, 9),
    ] {
        let a = pseudo_random_matrix(m, n, seed);
        assert_valid_thin_svd(&a, m, n, 1e-10);
    }
}

#[test]
fn svd_small_matches_closed_form_singular_values() {
    // A = [[4, 0], [3, -5]]: A^T A = [[25, -15], [-15, 25]] has eigenvalues 40 and 10.
    let a = [4.0, 0.0, 3.0, -5.0];
    let (_, s, _) = svd_small(&a, 2, 2);
    assert!((s[0] - 40.0_f64.sqrt()).abs() < 1e-13, "s0 = {}", s[0]);
    assert!((s[1] - 10.0_f64.sqrt()).abs() < 1e-13, "s1 = {}", s[1]);
    // The transpose has the same singular values.
    let a_t = [4.0, 3.0, 0.0, -5.0];
    let (_, s_t, _) = svd_small(&a_t, 2, 2);
    assert!((s_t[0] - s[0]).abs() < 1e-13 && (s_t[1] - s[1]).abs() < 1e-13);
    // A diagonal matrix: singular values are the sorted absolute diagonal.
    let diag = [3.0, 0.0, 0.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.5];
    let (_, s_diag, _) = svd_small(&diag, 3, 3);
    assert_eq!(s_diag.len(), 3);
    assert!((s_diag[0] - 4.0).abs() < 1e-14);
    assert!((s_diag[1] - 3.0).abs() < 1e-14);
    assert!((s_diag[2] - 0.5).abs() < 1e-14);
}

#[test]
fn svd_small_of_an_all_ones_matrix_has_one_singular_value() {
    for &(m, n) in &[(4usize, 6usize), (6, 4), (3, 3)] {
        let a = vec![1.0; m * n];
        let (u, s, vt) = svd_small(&a, m, n);
        let expected = ((m * n) as f64).sqrt();
        assert!((s[0] - expected).abs() < 1e-13, "{m}x{n}: s0 = {}", s[0]);
        assert!(s[1..].iter().all(|&v| v < 1e-13), "{m}x{n}: {s:?}");
        let k = m.min(n);
        // The zero singular values must not leave zero or non-orthogonal vectors.
        assert!(column_orthonormality_error(&u, m, k) < 1e-12, "{m}x{n}: U");
        assert!(row_orthonormality_error(&vt, k, n) < 1e-12, "{m}x{n}: Vt");
        assert!(svd_reconstruction_error(&a, (m, n, k), &u, &s, &vt) < 1e-12);
    }
}

#[test]
fn svd_small_of_a_rank_deficient_matrix_keeps_orthonormal_factors() {
    // Rank 2 by construction: third row = first row + second row.
    let a = [
        1.0, 2.0, 0.5, -1.0, //
        0.0, 1.0, 2.0, 3.0, //
        1.0, 3.0, 2.5, 2.0,
    ];
    let (u, s, vt) = svd_small(&a, 3, 4);
    assert!(s[0] > s[1] && s[1] > 1e-3 && s[2] < 1e-13, "{s:?}");
    assert!(column_orthonormality_error(&u, 3, 3) < 1e-12);
    assert!(row_orthonormality_error(&vt, 3, 4) < 1e-12);
    assert!(svd_reconstruction_error(&a, (3, 4, 3), &u, &s, &vt) < 1e-12);
}

#[test]
fn svd_small_zero_and_degenerate_inputs() {
    let (u, s, vt) = svd_small(&[0.0; 12], 3, 4);
    assert!(s.iter().all(|&v| v == 0.0));
    assert!(column_orthonormality_error(&u, 3, 3) < 1e-14);
    assert!(row_orthonormality_error(&vt, 3, 4) < 1e-14);
    let (u, s, vt) = svd_small(&[], 0, 4);
    assert!(u.is_empty() && s.is_empty() && vt.is_empty());
    let (u, s, vt) = svd_small(&[], 4, 0);
    assert!(u.is_empty() && s.is_empty() && vt.is_empty());
    let (u, s, vt) = svd_small(&[-3.0], 1, 1);
    assert_eq!(s, vec![3.0]);
    assert!((u[0] * s[0] * vt[0] + 3.0).abs() < 1e-15);
}

#[test]
fn svd_small_is_scale_robust() {
    // Squared column norms would overflow / underflow without internal scaling.
    for &scale in &[1e150_f64, 1e-150, 1e200] {
        let base = pseudo_random_matrix(3, 4, 21);
        let a: Vec<f64> = base.iter().map(|v| v * scale).collect();
        let (u, s, vt) = svd_small(&a, 3, 4);
        let (_, s_base, _) = svd_small(&base, 3, 4);
        for (scaled, unit) in s.iter().zip(s_base.iter()) {
            assert!(
                (scaled / scale - unit).abs() < 1e-13,
                "scale {scale}: {scaled} vs {unit}"
            );
        }
        let rescaled_error = svd_reconstruction_error(&a, (3, 4, 3), &u, &s, &vt) / scale;
        assert!(rescaled_error < 1e-13, "scale {scale}: {rescaled_error}");
    }
}

#[test]
fn qr_thin_returns_orthonormal_columns_for_a_dependent_input() {
    // Third column = first + second: the old routine normalised rounding noise
    // and returned a column that is not orthogonal to the first two.
    let a = [
        1.0, 2.0, 3.0, //
        2.0, -1.0, 1.0, //
        0.5, 0.0, 0.5, //
        -1.0, 1.0, 0.0, //
        3.0, 0.5, 3.5,
    ];
    let q = qr_thin(&a, 5, 3);
    assert!(column_orthonormality_error(&q, 5, 3) < 1e-12);
    // The first two columns still span the input's column space.
    for col in 0..3 {
        let mut residual = 0.0_f64;
        for i in 0..5 {
            let projected: f64 = (0..2)
                .map(|p| {
                    let coeff: f64 = (0..5).map(|r| q[r * 3 + p] * a[r * 3 + col]).sum();
                    coeff * q[i * 3 + p]
                })
                .sum();
            residual = residual.max((projected - a[i * 3 + col]).abs());
        }
        assert!(residual < 1e-12, "column {col} left the span: {residual}");
    }
}

#[test]
fn randomized_svd_of_an_all_ones_matrix_has_singular_value_sqrt_mn() {
    // The old code returned sqrt(m) because the k x n sketch was never rotated.
    let (m, n) = (6usize, 8usize);
    let a = vec![1.0; m * n];
    let res = LowRankApproximation::randomized_svd(&a, m, n, 1, 2);
    assert_eq!(res.u.len(), m);
    assert_eq!(res.vt.len(), n);
    assert!(
        (res.s[0] - ((m * n) as f64).sqrt()).abs() < 1e-12,
        "s0 = {}",
        res.s[0]
    );
    assert!(svd_reconstruction_error(&a, (m, n, 1), &res.u, &res.s, &res.vt) < 1e-12);
}

#[test]
fn randomized_svd_recovers_a_rank_two_matrix_with_a_larger_sketch() {
    // k = 3 exceeds the rank 2, so the sketch is rank deficient.
    let (m, n) = (7usize, 5usize);
    let left = pseudo_random_matrix(m, 2, 31);
    let right = pseudo_random_matrix(2, n, 32);
    let mut a = vec![0.0; m * n];
    for i in 0..m {
        for j in 0..n {
            a[i * n + j] = (0..2).map(|l| left[i * 2 + l] * right[l * n + j]).sum();
        }
    }
    let res = LowRankApproximation::randomized_svd(&a, m, n, 3, 2);
    assert_eq!((res.k, res.m, res.n), (3, m, n));
    assert_eq!(res.u.len(), m * 3);
    assert_eq!(res.vt.len(), 3 * n);
    let (_, exact, _) = svd_small(&a, m, n);
    assert!(
        (res.s[0] - exact[0]).abs() < 1e-10,
        "{} vs {}",
        res.s[0],
        exact[0]
    );
    assert!(
        (res.s[1] - exact[1]).abs() < 1e-10,
        "{} vs {}",
        res.s[1],
        exact[1]
    );
    assert!(res.s[2] < 1e-10, "third singular value {}", res.s[2]);
    assert!(svd_reconstruction_error(&a, (m, n, 3), &res.u, &res.s, &res.vt) < 1e-10);
    assert!(column_orthonormality_error(&res.u, m, 3) < 1e-10);
}

#[test]
fn nystrom_reconstructs_a_low_rank_psd_matrix() {
    // A = B B^T with B 6 x 2 has rank 2; three sampled columns span its range.
    let n = 6usize;
    let b = pseudo_random_matrix(n, 2, 41);
    let mut a = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            a[i * n + j] = (0..2).map(|l| b[i * 2 + l] * b[j * 2 + l]).sum();
        }
    }
    let res = LowRankApproximation::nystrom(&a, n, 3);
    assert_eq!((res.m, res.n, res.k), (n, n, 3));
    assert_eq!(res.u.len(), n * 3);
    assert_eq!(res.vt.len(), 3 * n);
    assert!(res.s.windows(2).all(|w| w[0] >= w[1]), "{:?}", res.s);
    assert!(column_orthonormality_error(&res.u, n, 3) < 1e-9);
    assert!(svd_reconstruction_error(&a, (n, n, 3), &res.u, &res.s, &res.vt) < 1e-9);
    // The eigenvalues of the approximation are the eigenvalues of A.
    let dense: Vec<Vec<f64>> = (0..n).map(|i| a[i * n..(i + 1) * n].to_vec()).collect();
    let (eigenvalues, _) = symmetric_eigen_n(&dense);
    assert!((res.s[0] - eigenvalues[n - 1]).abs() < 1e-9);
    assert!((res.s[1] - eigenvalues[n - 2]).abs() < 1e-9);
    assert!(res.s[2].abs() < 1e-9);
    // Degenerate request.
    let empty = LowRankApproximation::nystrom(&a, n, 0);
    assert_eq!(empty.k, 0);
    assert!(empty.u.is_empty() && empty.s.is_empty() && empty.vt.is_empty());
}

/// Contract the factors of a Tucker decomposition back into the full tensor.
fn tucker_to_full(res: &TensorDecompResult, shape: &[usize]) -> Vec<f64> {
    let ranks: Vec<usize> = res
        .factors
        .iter()
        .zip(shape.iter())
        .map(|(f, &dim)| f.len() / dim)
        .collect();
    let (di, dj, dk) = (shape[0], shape[1], shape[2]);
    let (r1, r2, r3) = (ranks[0], ranks[1], ranks[2]);
    assert_eq!(res.core.len(), r1 * r2 * r3);
    let mut full = vec![0.0; di * dj * dk];
    for i in 0..di {
        for j in 0..dj {
            for k in 0..dk {
                let mut acc = 0.0;
                for p in 0..r1 {
                    for q in 0..r2 {
                        for s in 0..r3 {
                            acc += res.core[(p * r2 + q) * r3 + s]
                                * res.factors[0][i * r1 + p]
                                * res.factors[1][j * r2 + q]
                                * res.factors[2][k * r3 + s];
                        }
                    }
                }
                full[(i * dj + j) * dk + k] = acc;
            }
        }
    }
    full
}

#[test]
fn tucker_hosvd_with_full_ranks_reconstructs_the_tensor() {
    for shape in [[4usize, 3, 2], [3, 3, 3], [6, 2, 2]] {
        let total: usize = shape.iter().product();
        let x = pseudo_random_matrix(total, 1, 51);
        let res = TensorDecomposition::tucker_hosvd(&x, &shape, &shape);
        let full = tucker_to_full(&res, &shape);
        let error = full
            .iter()
            .zip(x.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        assert!(
            error < 1e-10,
            "shape {shape:?}: reconstruction error {error}"
        );
    }
}

#[test]
fn tucker_hosvd_factors_are_orthonormal_and_ranks_are_clamped() {
    // A mode-3 unfolding of a 4 x 3 x 2 tensor has only 2 left singular
    // vectors, and the mode-1 unfolding of a 6 x 2 x 2 tensor only 4.
    let shape = [4usize, 3, 2];
    let x = pseudo_random_matrix(24, 1, 52);
    let res = TensorDecomposition::tucker_hosvd(&x, &shape, &[10, 10, 10]);
    assert_eq!(res.factors[0].len(), 4 * 4);
    assert_eq!(res.factors[1].len(), 3 * 3);
    assert_eq!(res.factors[2].len(), 2 * 2);
    assert_eq!(res.core.len(), 4 * 3 * 2);
    for (factor, &(rows, cols)) in res.factors.iter().zip(&[(4usize, 4usize), (3, 3), (2, 2)]) {
        assert!(column_orthonormality_error(factor, rows, cols) < 1e-12);
    }
    let tall =
        TensorDecomposition::tucker_hosvd(&pseudo_random_matrix(24, 1, 53), &[6, 2, 2], &[6, 2, 2]);
    assert_eq!(tall.factors[0].len(), 6 * 4);
    assert!(column_orthonormality_error(&tall.factors[0], 6, 4) < 1e-12);
}

#[test]
fn tucker_hosvd_captures_a_rank_one_tensor_with_unit_ranks() {
    let (a, b, c) = ([1.0, -2.0, 0.5, 3.0], [2.0, 1.0, -1.0], [1.5, -0.5]);
    let mut x = Vec::new();
    for ai in a {
        for bj in b {
            for ck in c {
                x.push(ai * bj * ck);
            }
        }
    }
    let shape = [4usize, 3, 2];
    let res = TensorDecomposition::tucker_hosvd(&x, &shape, &[1, 1, 1]);
    assert_eq!(res.core.len(), 1);
    let full = tucker_to_full(&res, &shape);
    for (got, want) in full.iter().zip(x.iter()) {
        assert!((got - want).abs() < 1e-12, "got {got}, want {want}");
    }
    // The mode-1 factor is the normalised first vector (up to sign).
    let norm_a = a.iter().map(|v| v * v).sum::<f64>().sqrt();
    let sign = if res.factors[0][0] * a[0] >= 0.0 {
        1.0
    } else {
        -1.0
    };
    for (u, v) in res.factors[0].iter().zip(a.iter()) {
        assert!((sign * u - v / norm_a).abs() < 1e-12);
    }
}

/// Contract TT cores (each `r_prev x n_k x r_k`, row-major) into the full tensor.
fn tt_to_full(cores: &[Vec<f64>], shape: &[usize]) -> Vec<f64> {
    let mut partial = vec![1.0_f64];
    let (mut r_prev, mut combos) = (1usize, 1usize);
    for (core, &n_k) in cores.iter().zip(shape.iter()) {
        let r_next = core.len() / (r_prev * n_k);
        let mut next = vec![0.0; combos * n_k * r_next];
        for c in 0..combos {
            for i in 0..n_k {
                for b in 0..r_next {
                    let mut acc = 0.0;
                    for a in 0..r_prev {
                        acc += partial[c * r_prev + a] * core[(a * n_k + i) * r_next + b];
                    }
                    next[(c * n_k + i) * r_next + b] = acc;
                }
            }
        }
        partial = next;
        combos *= n_k;
        r_prev = r_next;
    }
    assert_eq!(r_prev, 1, "the last TT rank must be 1");
    partial
}

#[test]
fn tensor_train_with_full_ranks_reconstructs_the_tensor() {
    for shape in [
        vec![3usize, 4, 2],
        vec![2, 3, 2, 3],
        vec![5],
        vec![2, 3],
        vec![4, 1, 3],
    ] {
        let total: usize = shape.iter().product();
        let x = pseudo_random_matrix(total, 1, 61);
        let res = TensorDecomposition::tensor_train(&x, &shape, 100, 0.0);
        assert_eq!(res.factors.len(), shape.len());
        let full = tt_to_full(&res.factors, &shape);
        assert_eq!(full.len(), total);
        let error = full
            .iter()
            .zip(x.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        assert!(
            error < 1e-10,
            "shape {shape:?}: reconstruction error {error}"
        );
    }
}

#[test]
fn tensor_train_of_a_rank_one_tensor_with_unit_ranks_is_exact() {
    // The last core has to carry the norm of the remainder.
    let (a, b, c) = ([1.0, -2.0, 0.5], [2.0, 1.0, -1.0, 4.0], [1.5, -0.5]);
    let mut x = Vec::new();
    for ai in a {
        for bj in b {
            for ck in c {
                x.push(ai * bj * ck);
            }
        }
    }
    let shape = [3usize, 4, 2];
    let res = TensorDecomposition::tensor_train(&x, &shape, 1, 0.0);
    assert_eq!(res.factors[0].len(), 3);
    assert_eq!(res.factors[1].len(), 4);
    assert_eq!(res.factors[2].len(), 2);
    let full = tt_to_full(&res.factors, &shape);
    for (got, want) in full.iter().zip(x.iter()) {
        assert!((got - want).abs() < 1e-12, "got {got}, want {want}");
    }
}

#[test]
fn tensor_train_max_rank_bounds_the_internal_ranks() {
    let shape = [3usize, 4, 3, 2];
    let total: usize = shape.iter().product();
    let x = pseudo_random_matrix(total, 1, 62);
    let res = TensorDecomposition::tensor_train(&x, &shape, 2, 0.0);
    let mut r_prev = 1usize;
    for (core, &n_k) in res.factors.iter().zip(shape.iter()) {
        let r_next = core.len() / (r_prev * n_k);
        assert!(r_next <= 2, "internal rank {r_next} exceeds max_rank");
        r_prev = r_next;
    }
    assert_eq!(r_prev, 1);
    assert_eq!(
        res.core.len(),
        res.factors.iter().map(Vec::len).sum::<usize>()
    );
}
