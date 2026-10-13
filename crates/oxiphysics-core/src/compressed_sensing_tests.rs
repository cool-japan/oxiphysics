// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Unit tests for the `compressed_sensing` module (sparse recovery, OMP,
//! K-SVD, RIP, MRI sampling).
//!
//! Split out of `compressed_sensing.rs` to keep that file under the 2000-line limit.
//! Included via `#[path]` as the body of `compressed_sensing::tests`.

use super::*;

// ---- soft_threshold ----

#[test]
fn test_soft_threshold_positive_above() {
    assert!((soft_threshold(3.0, 1.0) - 2.0).abs() < 1e-12);
}

#[test]
fn test_soft_threshold_positive_below() {
    assert_eq!(soft_threshold(0.5, 1.0), 0.0);
}

#[test]
fn test_soft_threshold_negative_above() {
    assert!((soft_threshold(-3.0, 1.0) + 2.0).abs() < 1e-12);
}

#[test]
fn test_soft_threshold_zero() {
    assert_eq!(soft_threshold(0.0, 1.0), 0.0);
}

#[test]
fn test_soft_threshold_exact_boundary() {
    assert_eq!(soft_threshold(1.0, 1.0), 0.0);
    assert_eq!(soft_threshold(-1.0, 1.0), 0.0);
}

#[test]
fn test_soft_threshold_zero_lambda() {
    assert!((soft_threshold(5.0, 0.0) - 5.0).abs() < 1e-12);
}

// ---- nyquist_rate ----

#[test]
fn test_nyquist_rate_basic() {
    assert!((nyquist_rate(1000.0) - 2000.0).abs() < 1e-9);
}

#[test]
fn test_nyquist_rate_zero() {
    assert_eq!(nyquist_rate(0.0), 0.0);
}

// ---- compression_ratio ----

#[test]
fn test_compression_ratio_half() {
    assert!((compression_ratio(100, 50) - 0.5).abs() < 1e-12);
}

#[test]
fn test_compression_ratio_zero_n() {
    assert_eq!(compression_ratio(0, 10), 0.0);
}

#[test]
fn test_compression_ratio_full() {
    assert!((compression_ratio(10, 10) - 1.0).abs() < 1e-12);
}

// ---- l2_norm / normalise ----

#[test]
fn test_l2_norm_known() {
    let x = vec![3.0, 4.0];
    assert!((l2_norm(&x) - 5.0).abs() < 1e-12);
}

#[test]
fn test_normalise_unit_vector() {
    let mut v = vec![3.0, 0.0, 4.0];
    normalise(&mut v);
    assert!((l2_norm(&v) - 1.0).abs() < 1e-12);
}

#[test]
fn test_normalise_zero_vector_unchanged() {
    let mut v = vec![0.0, 0.0, 0.0];
    normalise(&mut v);
    assert_eq!(v, vec![0.0, 0.0, 0.0]);
}

// ---- mat_vec / mat_transpose_vec ----

#[test]
fn test_mat_vec_identity() {
    let a = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
    let x = vec![3.0, 7.0];
    let y = mat_vec(&a, &x);
    assert!((y[0] - 3.0).abs() < 1e-12);
    assert!((y[1] - 7.0).abs() < 1e-12);
}

#[test]
fn test_mat_transpose_vec_basic() {
    let a = vec![vec![1.0, 2.0, 3.0]]; // 1×3
    let x = vec![2.0]; // length 1
    let y = mat_transpose_vec(&a, &x);
    assert_eq!(y.len(), 3);
    assert!((y[0] - 2.0).abs() < 1e-12);
    assert!((y[1] - 4.0).abs() < 1e-12);
    assert!((y[2] - 6.0).abs() < 1e-12);
}

// ---- DctBasis ----

#[test]
fn test_dct_roundtrip() {
    let n = 8;
    let basis = DctBasis::new(n);
    let signal: Vec<f64> = (0..n).map(|i| (i as f64).sin()).collect();
    let coeffs = basis.transform(&signal);
    let recovered = basis.inverse(&coeffs);
    for (a, b) in signal.iter().zip(recovered.iter()) {
        assert!((a - b).abs() < 1e-10, "DCT roundtrip mismatch: {a} vs {b}");
    }
}

#[test]
fn test_dct_dc_component() {
    // Constant signal → only DC (k=0) is non-zero.
    let n = 8;
    let basis = DctBasis::new(n);
    let signal = vec![1.0_f64; n];
    let coeffs = basis.transform(&signal);
    for (k, &c) in coeffs.iter().enumerate().skip(1) {
        assert!(
            c.abs() < 1e-10,
            "non-DC coefficient k={k} should be ~0, got {}",
            c
        );
    }
    assert!(coeffs[0].abs() > 0.5, "DC component should be non-zero");
}

#[test]
fn test_dct_length_preserved() {
    let n = 16;
    let basis = DctBasis::new(n);
    let signal: Vec<f64> = vec![1.0; n];
    let coeffs = basis.transform(&signal);
    assert_eq!(coeffs.len(), n);
}

#[test]
fn test_dct_energy_preservation() {
    // Parseval's theorem: ||x||^2 ≈ ||X||^2 for orthonormal DCT.
    let n = 8;
    let basis = DctBasis::new(n);
    let signal: Vec<f64> = (0..n).map(|i| i as f64).collect();
    let coeffs = basis.transform(&signal);
    let e_signal: f64 = signal.iter().map(|x| x * x).sum();
    let e_coeffs: f64 = coeffs.iter().map(|x| x * x).sum();
    assert!((e_signal - e_coeffs).abs() / (e_signal + 1.0) < 1e-10);
}

#[test]
fn test_dct_new_n() {
    let basis = DctBasis::new(4);
    assert_eq!(basis.n, 4);
}

#[test]
fn test_dct_truncate_keeps_k_largest() {
    let basis = DctBasis::new(8);
    let coeffs = vec![1.0, 5.0, 0.1, 3.0, 0.0, 2.0, 0.0, 0.0];
    let truncated = basis.truncate(&coeffs, 2);
    let nonzero = truncated.iter().filter(|&&v| v.abs() > 1e-14).count();
    assert_eq!(nonzero, 2, "truncate(k=2) should leave 2 non-zeros");
    assert!((truncated[1] - 5.0).abs() < 1e-12, "5.0 should be kept");
    assert!((truncated[3] - 3.0).abs() < 1e-12, "3.0 should be kept");
}

// ---- RandomMeasurementMatrix ----

#[test]
fn test_random_measurement_matrix_dimensions() {
    let mat = RandomMeasurementMatrix::generate_gaussian(10, 20);
    assert_eq!(mat.m, 10);
    assert_eq!(mat.n, 20);
    assert_eq!(mat.matrix.len(), 10);
    assert_eq!(mat.matrix[0].len(), 20);
}

#[test]
fn test_measurement_output_length() {
    let mat = RandomMeasurementMatrix::generate_gaussian(5, 10);
    let x = vec![1.0_f64; 10];
    let y = mat.measure(&x);
    assert_eq!(y.len(), 5);
}

#[test]
fn test_measurement_linearity() {
    let mat = RandomMeasurementMatrix::generate_gaussian(5, 8);
    let x1: Vec<f64> = (0..8).map(|i| i as f64).collect();
    let x2: Vec<f64> = (0..8).map(|i| (8 - i) as f64).collect();
    let y1 = mat.measure(&x1);
    let y2 = mat.measure(&x2);
    let y_sum: Vec<f64> = x1.iter().zip(x2.iter()).map(|(a, b)| a + b).collect();
    let y_direct = mat.measure(&y_sum);
    for (a, b) in y_direct
        .iter()
        .zip(y1.iter().zip(y2.iter()).map(|(a, b)| a + b))
    {
        assert!((a - b).abs() < 1e-10, "linearity: {a} vs {b}");
    }
}

#[test]
fn test_bernoulli_matrix_dimensions() {
    let mat = RandomMeasurementMatrix::generate_bernoulli(8, 16);
    assert_eq!(mat.m, 8);
    assert_eq!(mat.n, 16);
}

#[test]
fn test_bernoulli_entries_are_plus_minus_scale() {
    let m = 5;
    let n = 10;
    let mat = RandomMeasurementMatrix::generate_bernoulli(m, n);
    let scale = 1.0 / (m as f64).sqrt();
    for row in &mat.matrix {
        for &v in row {
            let diff = (v.abs() - scale).abs();
            assert!(diff < 1e-12, "Bernoulli entry |{v}| ≠ {scale}");
        }
    }
}

// ---- BasisPursuit (ISTA) ----

#[test]
fn test_ista_trivial_identity() {
    // A = I_4, b = [1,0,0,0], lambda small → recover e_0
    let a: Vec<Vec<f64>> = (0..4)
        .map(|i| (0..4).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
        .collect();
    let b = vec![1.0, 0.0, 0.0, 0.0];
    let x = BasisPursuit::solve_lasso(&a, &b, 1e-4, 200);
    assert_eq!(x.len(), 4);
    assert!((x[0] - 1.0).abs() < 0.05, "x[0] should be ~1, got {}", x[0]);
    assert!(x[1].abs() < 0.05);
}

#[test]
fn test_ista_empty_input() {
    let x = BasisPursuit::solve_lasso(&[], &[], 1.0, 100);
    assert!(x.is_empty());
}

#[test]
fn test_ista_sparse_recovery() {
    // A = 4×4 identity, sparse signal [3, 0, 0, 0]
    let a: Vec<Vec<f64>> = (0..4)
        .map(|i| (0..4).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
        .collect();
    let b = vec![3.0, 0.0, 0.0, 0.0];
    let x = BasisPursuit::solve_lasso(&a, &b, 0.01, 300);
    assert!((x[0] - 3.0).abs() < 0.1, "x[0] ≈ 3, got {}", x[0]);
}

// ---- FISTA ----

#[test]
fn test_fista_identity_recovery() {
    let a: Vec<Vec<f64>> = (0..4)
        .map(|i| (0..4).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
        .collect();
    let b = vec![0.0, 2.5, 0.0, 0.0];
    let x = BasisPursuit::solve_fista(&a, &b, 1e-3, 300);
    assert!((x[1] - 2.5).abs() < 0.05, "x[1] ≈ 2.5, got {}", x[1]);
}

#[test]
fn test_fista_empty_input() {
    let x = BasisPursuit::solve_fista(&[], &[], 1.0, 100);
    assert!(x.is_empty());
}

#[test]
fn test_fista_objective_decreases() {
    let a: Vec<Vec<f64>> = (0..3)
        .map(|i| (0..3).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
        .collect();
    let b = vec![1.0, -1.0, 0.5];
    let lambda = 0.1;
    let x0 = vec![0.0_f64; 3];
    let obj0 = BasisPursuit::objective(&a, &b, &x0, lambda);
    let x_hat = BasisPursuit::solve_fista(&a, &b, lambda, 100);
    let obj1 = BasisPursuit::objective(&a, &b, &x_hat, lambda);
    assert!(
        obj1 <= obj0 + 1e-10,
        "FISTA should decrease objective: {obj1} > {obj0}"
    );
}

// ---- OMP ----

#[test]
fn test_omp_exact_1_sparse() {
    // Identity system with 1-sparse signal
    let a: Vec<Vec<f64>> = (0..4)
        .map(|i| (0..4).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
        .collect();
    let b = vec![0.0, 5.0, 0.0, 0.0];
    let omp = OrthogonalMatchingPursuit::new(1);
    let x = omp.solve(&a, &b);
    assert_eq!(x.len(), 4);
    assert!((x[1] - 5.0).abs() < 1e-10, "x[1] should be 5, got {}", x[1]);
}

#[test]
fn test_omp_exact_2_sparse() {
    let a: Vec<Vec<f64>> = (0..4)
        .map(|i| (0..4).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
        .collect();
    let b = vec![2.0, 0.0, 7.0, 0.0];
    let omp = OrthogonalMatchingPursuit::new(2);
    let x = omp.solve(&a, &b);
    assert!((x[0] - 2.0).abs() < 1e-8);
    assert!((x[2] - 7.0).abs() < 1e-8);
}

#[test]
fn test_omp_empty_input() {
    let omp = OrthogonalMatchingPursuit::new(3);
    let x = omp.solve(&[], &[]);
    assert!(x.is_empty());
}

#[test]
fn test_omp_new() {
    let omp = OrthogonalMatchingPursuit::new(5);
    assert_eq!(omp.max_k, 5);
}

#[test]
fn test_omp_residual_decreases() {
    // Random orthogonal system: with enough sparsity the residual should shrink.
    let a: Vec<Vec<f64>> = (0..6)
        .map(|i| (0..6).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
        .collect();
    let b = vec![1.0, -1.0, 2.0, 0.0, -2.0, 0.0];
    let omp = OrthogonalMatchingPursuit::new(4);
    let x = omp.solve(&a, &b);
    // Residual after recovery
    let residual: f64 = b
        .iter()
        .enumerate()
        .map(|(i, &bi)| {
            let ax: f64 = a[i].iter().zip(x.iter()).map(|(aij, xj)| aij * xj).sum();
            (bi - ax).powi(2)
        })
        .sum::<f64>()
        .sqrt();
    assert!(residual < 1e-8, "residual should be tiny, got {residual}");
}

#[test]
fn test_omp_support_length() {
    let a: Vec<Vec<f64>> = (0..5)
        .map(|i| (0..5).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
        .collect();
    let b = vec![1.0, 0.0, 3.0, 0.0, 2.0];
    let omp = OrthogonalMatchingPursuit::new(3);
    let supp = omp.support(&a, &b);
    assert_eq!(
        supp.len(),
        3,
        "support should have 3 elements, got {}",
        supp.len()
    );
}

// ---- SparsityMetrics ----

#[test]
fn test_l0_norm_basic() {
    let x = vec![0.0, 1.0, 0.0, -2.0, 0.001];
    assert_eq!(SparsityMetrics::l0_norm(&x, 0.5), 2); // 1.0 and -2.0
}

#[test]
fn test_l0_norm_all_zero() {
    let x = vec![0.0, 0.0, 0.0];
    assert_eq!(SparsityMetrics::l0_norm(&x, 1e-6), 0);
}

#[test]
fn test_l1_norm_basic() {
    let x = vec![1.0, -2.0, 3.0];
    assert!((SparsityMetrics::l1_norm(&x) - 6.0).abs() < 1e-12);
}

#[test]
fn test_l1_norm_empty() {
    assert_eq!(SparsityMetrics::l1_norm(&[]), 0.0);
}

#[test]
fn test_l2_norm_sparsity() {
    let x = vec![3.0, 4.0];
    assert!((SparsityMetrics::l2_norm(&x) - 5.0).abs() < 1e-12);
}

#[test]
fn test_gini_sparse_is_near_one() {
    let x = vec![0.0, 0.0, 0.0, 10.0]; // very sparse
    let g = SparsityMetrics::gini(&x);
    assert!(g > 0.6, "sparse signal should have high Gini, got {g}");
}

#[test]
fn test_gini_uniform_is_near_zero() {
    let x = vec![1.0, 1.0, 1.0, 1.0]; // flat = dense
    let g = SparsityMetrics::gini(&x);
    assert!(g < 0.1, "uniform signal should have low Gini, got {g}");
}

#[test]
fn test_coherence_identity_is_zero() {
    // Columns of the identity are orthogonal → coherence = 0
    let a: Vec<Vec<f64>> = (0..4)
        .map(|i| (0..4).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
        .collect();
    let mu = SparsityMetrics::coherence(&a);
    assert!(mu < 1e-12, "identity coherence should be 0, got {mu}");
}

#[test]
fn test_coherence_empty() {
    assert_eq!(SparsityMetrics::coherence(&[]), 0.0);
}

#[test]
fn test_coherence_collinear_columns() {
    // Two identical columns → coherence = 1
    let a = vec![vec![1.0, 1.0], vec![0.0, 0.0]];
    let mu = SparsityMetrics::coherence(&a);
    assert!(
        (mu - 1.0).abs() < 1e-10,
        "collinear columns → coherence=1, got {mu}"
    );
}

#[test]
fn test_babel_function_identity() {
    let a: Vec<Vec<f64>> = (0..4)
        .map(|i| (0..4).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
        .collect();
    let babel = SparsityMetrics::babel_function(&a, 2);
    assert!(babel < 1e-12, "identity babel should be 0, got {babel}");
}

// ---- RecoveryGuarantee ----

#[test]
fn test_exact_recovery_condition_sufficient_measurements() {
    // k=1, n=100 → need m >= 2*ln(100) ≈ 9.2 → m=20 is sufficient
    assert!(RecoveryGuarantee::exact_recovery_condition(1, 20, 100));
}

#[test]
fn test_exact_recovery_condition_insufficient() {
    // k=50, n=100 → need m >= 100*ln(2) ≈ 69 → m=5 is not sufficient
    assert!(!RecoveryGuarantee::exact_recovery_condition(50, 5, 100));
}

#[test]
fn test_exact_recovery_condition_k_zero() {
    assert!(RecoveryGuarantee::exact_recovery_condition(0, 0, 100));
}

#[test]
fn test_rip_constant_identity() {
    // For the identity matrix (m=n), RIP constant should be ~0
    let n = 4;
    let a: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
        .collect();
    let delta = RecoveryGuarantee::rip_constant(&a, 1);
    assert!(
        delta < 1e-10,
        "identity RIP constant should be ~0, got {delta}"
    );
}

#[test]
fn test_rip_constant_empty() {
    assert_eq!(RecoveryGuarantee::rip_constant(&[], 1), 0.0);
}

#[test]
fn test_rip_measurement_lower_bound_nonzero() {
    let lb = RecoveryGuarantee::rip_measurement_lower_bound(5, 100);
    assert!(lb > 0, "lower bound should be positive");
}

#[test]
fn test_lasso_error_bound_positive() {
    let bound = RecoveryGuarantee::lasso_error_bound(0.1, 4, 0.1);
    assert!(bound > 0.0, "LASSO error bound should be positive");
}

// ---- gauss_solve ----

#[test]
fn test_gauss_solve_2x2() {
    // [2 1; 1 3] x = [5; 10] → x = [1, 3]
    let a = vec![vec![2.0, 1.0], vec![1.0, 3.0]];
    let b = vec![5.0, 10.0];
    let x = gauss_solve(&a, &b);
    assert!((x[0] - 1.0).abs() < 1e-10, "x[0]={}", x[0]);
    assert!((x[1] - 3.0).abs() < 1e-10, "x[1]={}", x[1]);
}

#[test]
fn test_gauss_solve_1x1() {
    let a = vec![vec![4.0]];
    let b = vec![8.0];
    let x = gauss_solve(&a, &b);
    assert!((x[0] - 2.0).abs() < 1e-12);
}

#[test]
fn test_gauss_solve_empty() {
    let x = gauss_solve(&[], &[]);
    assert!(x.is_empty());
}

// ---- MriCompressedSensing ----

#[test]
fn test_mri_cs_new() {
    let mri = MriCompressedSensing::new(64, 20);
    assert_eq!(mri.n, 64);
    assert_eq!(mri.m, 20);
}

#[test]
fn test_mri_kspace_indices_length() {
    let mri = MriCompressedSensing::new(32, 10);
    let idx = mri.sample_kspace_indices();
    assert_eq!(idx.len(), 10);
}

#[test]
fn test_mri_measurement_matrix_shape() {
    let mri = MriCompressedSensing::new(16, 8);
    let idx: Vec<usize> = (0..8).collect();
    let a = mri.build_measurement_matrix(&idx);
    assert_eq!(a.len(), 8);
    assert_eq!(a[0].len(), 16);
}

#[test]
fn test_psnr_identical_signals() {
    let s = vec![1.0, 2.0, 3.0];
    let psnr = MriCompressedSensing::psnr(&s, &s, 3.0);
    assert!(psnr.is_infinite(), "identical signals → PSNR = ∞");
}

#[test]
fn test_psnr_known_value() {
    let original = vec![1.0, 0.0];
    let reconstructed = vec![0.0, 0.0];
    let psnr = MriCompressedSensing::psnr(&original, &reconstructed, 1.0);
    assert!(psnr.is_finite(), "PSNR should be finite for non-identical");
}

// ---- SparseSignal ----

#[test]
fn test_sparse_signal_generate_sparsity() {
    let sig = SparseSignal::generate(20, 3, 1.0);
    assert_eq!(sig.len(), 20);
    let nnz = sig.iter().filter(|&&v| v.abs() > 1e-14).count();
    assert_eq!(
        nnz, 3,
        "generated signal should have exactly 3 non-zeros, got {nnz}"
    );
}

#[test]
fn test_sparse_signal_generate_length() {
    let sig = SparseSignal::generate(100, 5, 2.0);
    assert_eq!(sig.len(), 100);
}

#[test]
fn test_sparse_signal_relative_error_zero() {
    let s = vec![1.0, 2.0, 3.0];
    let err = SparseSignal::relative_error(&s, &s);
    assert!(
        err < 1e-12,
        "identical signals should have 0 relative error"
    );
}

#[test]
fn test_sparse_signal_support_error_identical() {
    let s = vec![0.0, 1.0, 0.0, 2.0];
    let err = SparseSignal::support_error(&s, &s, 0.5);
    assert_eq!(err, 0.0, "identical support → error = 0");
}

// ---- KSvd ----

#[test]
fn test_ksvd_new() {
    let k = KSvd::new(8, 2, 5);
    assert_eq!(k.n_atoms, 8);
    assert_eq!(k.sparsity, 2);
    assert_eq!(k.n_iter, 5);
}

#[test]
fn test_ksvd_reconstruct_zero_code() {
    let dict: Vec<Vec<f64>> = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
    let code = vec![0.0, 0.0];
    let rec = KSvd::reconstruct(&dict, &code);
    assert_eq!(rec, vec![0.0, 0.0]);
}

#[test]
fn test_ksvd_reconstruct_unit_code() {
    let dict: Vec<Vec<f64>> = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
    let code = vec![3.0, 5.0];
    let rec = KSvd::reconstruct(&dict, &code);
    assert!((rec[0] - 3.0).abs() < 1e-12);
    assert!((rec[1] - 5.0).abs() < 1e-12);
}

#[test]
fn test_ksvd_fit_returns_correct_shape() {
    // 4 signals of length 6, 4 atoms
    let signals: Vec<Vec<f64>> = (0..4)
        .map(|i| (0..6).map(|j| if j == i { 1.0 } else { 0.0 }).collect())
        .collect();
    let ksvd = KSvd::new(4, 1, 2);
    let dict = ksvd.fit(&signals);
    assert_eq!(dict.len(), 4, "dict should have 4 atoms");
    assert_eq!(dict[0].len(), 6, "each atom should have length 6");
}

/// Random unit-norm dictionary of `n_atoms` atoms in ℝᵈ and `n_signals`
/// signals, each a combination of `active` distinct atoms (seeded).
fn planted_problem(
    d: usize,
    n_atoms: usize,
    active: usize,
    n_signals: usize,
    seed: u64,
) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let mut rng = crate::StatRng::new(seed);
    let atoms: Vec<Vec<f64>> = (0..n_atoms)
        .map(|_| {
            let mut a: Vec<f64> = (0..d).map(|_| rng.next_normal()).collect();
            normalise(&mut a);
            a
        })
        .collect();
    let signals: Vec<Vec<f64>> = (0..n_signals)
        .map(|_| {
            let mut y = vec![0.0_f64; d];
            let mut chosen: Vec<usize> = Vec::new();
            while chosen.len() < active {
                let idx = ((rng.next_f64() * n_atoms as f64) as usize).min(n_atoms - 1);
                if !chosen.contains(&idx) {
                    chosen.push(idx);
                }
            }
            for &idx in &chosen {
                let c = rng.next_normal();
                for (yi, &ai) in y.iter_mut().zip(atoms[idx].iter()) {
                    *yi += c * ai;
                }
            }
            y
        })
        .collect();
    (atoms, signals)
}

#[test]
fn test_ksvd_default_seed_and_with_seed() {
    let k = KSvd::new(8, 2, 5);
    assert_eq!(k.seed, KSVD_DEFAULT_SEED);
    let k = k.with_seed(99);
    assert_eq!(k.seed, 99);
    assert_eq!((k.n_atoms, k.sparsity, k.n_iter), (8, 2, 5));
}

#[test]
fn test_ksvd_fit_is_deterministic_for_a_seed() {
    let (_, signals) = planted_problem(8, 12, 2, 80, 11);
    let a = KSvd::new(12, 2, 4).with_seed(5).fit(&signals);
    let b = KSvd::new(12, 2, 4).with_seed(5).fit(&signals);
    assert_eq!(a, b, "same seed must give an identical dictionary");
    let c = KSvd::new(12, 2, 4).with_seed(6).fit(&signals);
    assert_ne!(a, c, "a different seed must change the initialisation");
}

#[test]
fn test_ksvd_overcomplete_dictionary_is_not_capped_at_signal_dimension() {
    // 2× overcomplete: 12 atoms for signals of dimension 6.
    let (_, signals) = planted_problem(6, 12, 2, 60, 3);
    let dict = KSvd::new(12, 2, 3).fit(&signals);
    assert_eq!(dict.len(), 12, "all requested atoms must be learned");
    for atom in &dict {
        assert_eq!(atom.len(), 6, "atoms have the signal dimension");
        assert!((l2_norm(atom) - 1.0).abs() < 1e-9, "unit-norm atoms");
    }
}

#[test]
fn test_ksvd_update_stage_never_increases_representation_error() {
    // Refreshing the coefficients together with each atom makes the update
    // stage a monotone (rank-1 least-squares) step.  Stale coefficients
    // break this guarantee.
    for seed in [1_u64, 2, 3, 4] {
        let (_, signals) = planted_problem(10, 20, 3, 300, 100 + seed);
        let (_, history) = KSvd::new(20, 3, 8)
            .with_seed(seed)
            .fit_with_history(&signals);
        assert_eq!(history.len(), 8);
        for (it, h) in history.iter().enumerate() {
            assert!(
                h.after_update <= h.after_coding * (1.0 + 1e-9) + 1e-12,
                "seed {seed}, iteration {it}: update raised the error {} -> {}",
                h.after_coding,
                h.after_update
            );
        }
    }
}

#[test]
fn test_ksvd_error_decreases_on_planted_overcomplete_dictionary() {
    // 2× overcomplete planted dictionary (20 atoms in ℝ¹⁰), 3-sparse codes.
    let (_, signals) = planted_problem(10, 20, 3, 300, 2024);
    let (dict, history) = KSvd::new(20, 3, 10).with_seed(7).fit_with_history(&signals);
    assert_eq!(dict.len(), 20);
    let mut previous = f64::INFINITY;
    for (it, h) in history.iter().enumerate() {
        assert!(
            h.after_update <= previous * (1.0 + 1e-9) + 1e-12,
            "iteration {it}: error rose from {previous} to {}",
            h.after_update
        );
        previous = h.after_update;
    }
    let signal_norm = signals
        .iter()
        .map(|y| l2_norm(y).powi(2))
        .sum::<f64>()
        .sqrt();
    let first = history.first().map(|h| h.after_coding).unwrap_or(f64::NAN);
    let last = history.last().map(|h| h.after_update).unwrap_or(f64::NAN);
    assert!(
        last < 0.8 * first,
        "dictionary learning must reduce the error: {first} -> {last} (‖Y‖ = {signal_norm})"
    );
}

#[test]
fn test_ksvd_fit_with_history_matches_fit() {
    let (_, signals) = planted_problem(8, 12, 2, 80, 21);
    let ksvd = KSvd::new(12, 2, 3).with_seed(8);
    let (dict, history) = ksvd.fit_with_history(&signals);
    assert_eq!(dict, ksvd.fit(&signals));
    assert_eq!(history.len(), 3);
}

#[test]
fn test_ksvd_recovers_axis_aligned_dictionary() {
    // Signals are ±unit multiples of the 4 canonical axes of ℝ⁴: a dictionary
    // of 4 atoms with 1-sparse codes can represent them exactly.
    let mut signals: Vec<Vec<f64>> = Vec::new();
    for rep in 0..6 {
        for axis in 0..4 {
            let mut y = vec![0.0; 4];
            y[axis] = if rep % 2 == 0 { 1.5 } else { -0.7 };
            signals.push(y);
        }
    }
    let (dict, history) = KSvd::new(4, 1, 6).with_seed(1).fit_with_history(&signals);
    let last = history.last().map(|h| h.after_update).unwrap_or(f64::NAN);
    assert!(
        last < 1e-9,
        "exact 1-sparse structure must be recovered, got {last}"
    );
    // Each canonical axis is (up to sign) one of the learned atoms.
    for axis in 0..4 {
        let best = dict.iter().map(|a| a[axis].abs()).fold(0.0_f64, f64::max);
        assert!((best - 1.0).abs() < 1e-9, "axis {axis} not learned: {best}");
    }
}

#[test]
fn test_ksvd_rejects_degenerate_training_sets() {
    let ksvd = KSvd::new(4, 1, 2);
    assert!(ksvd.fit(&[]).is_empty(), "no signals");
    assert!(
        ksvd.fit(&[Vec::new(), Vec::new()]).is_empty(),
        "empty signals"
    );
    let ragged = vec![vec![1.0, 2.0, 3.0], vec![1.0, 2.0]];
    assert!(ksvd.fit(&ragged).is_empty(), "signals of different lengths");
    let (dict, history) = ksvd.fit_with_history(&[]);
    assert!(dict.is_empty() && history.is_empty());
}

#[test]
fn test_ksvd_encode_length() {
    let dict: Vec<Vec<f64>> = (0..4)
        .map(|i| (0..4).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
        .collect();
    let ksvd = KSvd::new(4, 1, 1);
    let signal = vec![0.0, 1.0, 0.0, 0.0];
    let code = ksvd.encode(&dict, &signal);
    assert_eq!(code.len(), 4, "code length should match n_atoms");
}

#[test]
fn test_spectral_norm_identity() {
    let a: Vec<Vec<f64>> = (0..4)
        .map(|i| (0..4).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
        .collect();
    let sn = spectral_norm(&a, 30);
    assert!(
        (sn - 1.0).abs() < 0.01,
        "spectral norm of I should be ~1, got {sn}"
    );
}

#[test]
fn test_spectral_norm_empty() {
    let sn = spectral_norm(&[], 10);
    assert_eq!(sn, 0.0);
}
