// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Unit tests of the eigenvalue solvers (child module of `eigenvalue`, so the
//! private numerical kernels are reachable through `super::*`).

use super::*;

// ── Utility helpers ───────────────────────────────────────────────────────

/// 2×2 symmetric positive-definite matrix: \[\[4,1\\],\[1,3\]]
fn mat2_spd() -> (Vec<f64>, usize) {
    (vec![4.0, 1.0, 1.0, 3.0], 2)
}

/// Identity matrix n×n
fn identity(n: usize) -> Vec<f64> {
    let mut a = vec![0.0_f64; n * n];
    for i in 0..n {
        a[i * n + i] = 1.0;
    }
    a
}

/// Diagonal matrix with given diagonal entries
fn diag_mat(diag: &[f64]) -> Vec<f64> {
    let n = diag.len();
    let mut a = vec![0.0_f64; n * n];
    for (i, &d) in diag.iter().enumerate() {
        a[i * n + i] = d;
    }
    a
}

// ── EigensolverConfig ─────────────────────────────────────────────────────

#[test]
fn test_eigensolver_config_default() {
    let cfg = EigensolverConfig::default();
    assert_eq!(cfg.max_iter, 500);
    assert_eq!(cfg.n_modes, 5);
    assert!(cfg.tolerance < 1e-9);
    assert_eq!(cfg.shift, 0.0);
}

#[test]
fn test_eigensolver_config_custom() {
    let cfg = EigensolverConfig {
        max_iter: 100,
        tolerance: 1e-6,
        shift: 2.0,
        n_modes: 3,
    };
    assert_eq!(cfg.n_modes, 3);
    assert_eq!(cfg.shift, 2.0);
}

// ── PowerIteration ────────────────────────────────────────────────────────

#[test]
fn test_power_iteration_identity() {
    let (a, n) = (identity(3), 3);
    let cfg = EigensolverConfig {
        max_iter: 200,
        tolerance: 1e-10,
        ..Default::default()
    };
    let pi = PowerIteration::new(a, n, cfg);
    let (lam, v) = pi.dominant();
    assert!(
        (lam - 1.0).abs() < 1e-6,
        "eigenvalue of I should be 1, got {lam}"
    );
    assert!((norm2(&v) - 1.0).abs() < 1e-10);
}

#[test]
fn test_power_iteration_diagonal() {
    let a = diag_mat(&[1.0, 5.0, 3.0]);
    let cfg = EigensolverConfig {
        max_iter: 1000,
        tolerance: 1e-10,
        ..Default::default()
    };
    let pi = PowerIteration::new(a, 3, cfg);
    let (lam, _) = pi.dominant();
    assert!(
        (lam - 5.0).abs() < 1e-5,
        "dominant eigenvalue should be 5, got {lam}"
    );
}

#[test]
fn test_power_iteration_2x2_spd() {
    // [[4,1],[1,3]] has eigenvalues ≈ 4.618, 2.382
    let (a, n) = mat2_spd();
    let cfg = EigensolverConfig {
        max_iter: 500,
        tolerance: 1e-10,
        ..Default::default()
    };
    let pi = PowerIteration::new(a, n, cfg);
    let (lam, v) = pi.dominant();
    assert!(
        lam > 4.0 && lam < 5.0,
        "dominant eigenvalue out of range: {lam}"
    );
    assert!((norm2(&v) - 1.0).abs() < 1e-10);
}

#[test]
fn test_power_iteration_deflation_count() {
    let a = diag_mat(&[1.0, 4.0, 9.0]);
    let cfg = EigensolverConfig {
        max_iter: 1000,
        tolerance: 1e-10,
        ..Default::default()
    };
    let pi = PowerIteration::new(a, 3, cfg);
    let results = pi.deflated_eigenvalues(2);
    assert_eq!(results.len(), 2);
}

#[test]
fn test_power_iteration_eigvec_unit_norm() {
    let a = diag_mat(&[2.0, 7.0, 3.0]);
    let cfg = EigensolverConfig {
        max_iter: 500,
        tolerance: 1e-10,
        ..Default::default()
    };
    let pi = PowerIteration::new(a, 3, cfg);
    let (_, v) = pi.dominant();
    assert!((norm2(&v) - 1.0).abs() < 1e-9);
}

// ── InverseIteration ──────────────────────────────────────────────────────

#[test]
fn test_inverse_iteration_identity() {
    let cfg = EigensolverConfig {
        max_iter: 300,
        tolerance: 1e-10,
        shift: 0.0,
        n_modes: 1,
    };
    let inv = InverseIteration::new(identity(3), 3, cfg);
    let (lam, _) = inv.solve().unwrap();
    assert!(
        (lam - 1.0).abs() < 1e-5,
        "smallest eigenvalue of I = 1, got {lam}"
    );
}

#[test]
fn test_inverse_iteration_diagonal() {
    let a = diag_mat(&[3.0, 1.0, 5.0]);
    let cfg = EigensolverConfig {
        max_iter: 500,
        tolerance: 1e-10,
        shift: 0.5,
        n_modes: 1,
    };
    let inv = InverseIteration::new(a, 3, cfg);
    let (lam, _) = inv.solve().unwrap();
    // Nearest eigenvalue to shift=0.5 is 1.0
    assert!((lam - 1.0).abs() < 1e-4, "expected 1.0, got {lam}");
}

#[test]
fn test_inverse_iteration_returns_unit_eigvec() {
    let a = diag_mat(&[2.0, 5.0, 8.0]);
    let cfg = EigensolverConfig {
        max_iter: 300,
        tolerance: 1e-10,
        shift: 1.5,
        n_modes: 1,
    };
    let inv = InverseIteration::new(a, 3, cfg);
    let (_, v) = inv.solve().unwrap();
    assert!((norm2(&v) - 1.0).abs() < 1e-9);
}

#[test]
fn test_inverse_iteration_2x2() {
    let (a, n) = mat2_spd();
    let cfg = EigensolverConfig {
        max_iter: 500,
        tolerance: 1e-10,
        shift: 2.0,
        n_modes: 1,
    };
    let inv = InverseIteration::new(a, n, cfg);
    let (lam, _) = inv.solve().unwrap();
    // Smallest eigenvalue ≈ 2.382
    assert!(lam > 2.0 && lam < 3.0, "expected ~2.38, got {lam}");
}

// ── LanczosMethod ─────────────────────────────────────────────────────────

#[test]
fn test_lanczos_identity_eigenvalues() {
    let cfg = EigensolverConfig {
        n_modes: 3,
        ..Default::default()
    };
    let lm = LanczosMethod::new(identity(4), 4, cfg);
    let eigs = lm.eigenvalues();
    for &e in &eigs {
        assert!((e - 1.0).abs() < 0.1, "identity eigenvalue {e} != 1");
    }
}

#[test]
fn test_lanczos_diagonal_sorted() {
    let a = diag_mat(&[5.0, 2.0, 8.0, 1.0]);
    let cfg = EigensolverConfig {
        n_modes: 4,
        max_iter: 300,
        tolerance: 1e-10,
        ..Default::default()
    };
    let lm = LanczosMethod::new(a, 4, cfg);
    let eigs = lm.eigenvalues();
    for w in eigs.windows(2) {
        assert!(w[0] <= w[1] + 1e-10, "eigenvalues not sorted");
    }
}

#[test]
fn test_lanczos_n_modes_limit() {
    let a = diag_mat(&[1.0, 2.0, 3.0, 4.0, 5.0]);
    let cfg = EigensolverConfig {
        n_modes: 3,
        ..Default::default()
    };
    let lm = LanczosMethod::new(a, 5, cfg);
    let eigs = lm.eigenvalues();
    assert!(eigs.len() <= 3);
}

#[test]
fn test_lanczos_2x2_spd() {
    let (a, n) = mat2_spd();
    let cfg = EigensolverConfig {
        n_modes: 2,
        max_iter: 300,
        tolerance: 1e-10,
        ..Default::default()
    };
    let lm = LanczosMethod::new(a, n, cfg);
    let eigs = lm.eigenvalues();
    assert!(!eigs.is_empty());
    // Both eigenvalues should be positive
    for &e in &eigs {
        assert!(e > 0.0, "eigenvalue should be positive, got {e}");
    }
}

// ── RayleighQuotient ──────────────────────────────────────────────────────

#[test]
fn test_rayleigh_quotient_identity() {
    let k = identity(3);
    let m = identity(3);
    let rq = RayleighQuotient::new(k, m, 3);
    let x = vec![1.0, 0.0, 0.0];
    let r = rq.evaluate(&x);
    assert!((r - 1.0).abs() < 1e-12);
}

#[test]
fn test_rayleigh_quotient_diagonal() {
    let k = diag_mat(&[1.0, 4.0, 9.0]);
    let m = identity(3);
    let rq = RayleighQuotient::new(k, m, 3);
    // x aligned with first mode
    let x = vec![1.0, 0.0, 0.0];
    let r = rq.evaluate(&x);
    assert!((r - 1.0).abs() < 1e-12);
}

#[test]
fn test_rayleigh_quotient_upper_bound() {
    let k = diag_mat(&[1.0, 4.0, 9.0]);
    let m = identity(3);
    let rq = RayleighQuotient::new(k, m, 3);
    // Any trial vector ρ(x) >= min eigenvalue = 1.0
    let x = vec![0.5, 0.5, 0.5f64.sqrt()];
    let r = rq.evaluate(&x);
    assert!(
        r >= 1.0 - 1e-10,
        "Rayleigh quotient must be >= min eigenvalue"
    );
}

#[test]
fn test_rayleigh_minimise() {
    let k = diag_mat(&[1.0, 4.0, 9.0]);
    let m = identity(3);
    let rq = RayleighQuotient::new(k, m, 3);
    let (rho_min, _) = rq.minimise(100);
    assert!(
        rho_min >= 1.0 - 1e-6,
        "minimised RQ must stay >= min eigenvalue"
    );
    assert!(rho_min <= 9.0 + 1e-6);
}

// ── GeneralizedEigenvalue ─────────────────────────────────────────────────

#[test]
fn test_generalized_eigenvalue_identity_mass() {
    // K = diag(1,4), M = I → eigenvalues = 1, 4
    let k = diag_mat(&[1.0, 4.0]);
    let m = identity(2);
    let cfg = EigensolverConfig {
        n_modes: 2,
        max_iter: 500,
        tolerance: 1e-10,
        ..Default::default()
    };
    let gev = GeneralizedEigenvalue::new(k, m, 2, cfg);
    let (eigs, _) = gev.solve_smallest().unwrap();
    assert_eq!(eigs.len(), 2);
    assert!((eigs[0] - 1.0).abs() < 1e-12, "expected 1, got {}", eigs[0]);
    assert!((eigs[1] - 4.0).abs() < 1e-12, "expected 4, got {}", eigs[1]);
}

#[test]
fn test_generalized_eigenvalue_cholesky_spd() {
    // 2×2 SPD mass matrix
    let m = vec![4.0, 2.0, 2.0, 3.0];
    let l = GeneralizedEigenvalue::cholesky(&m, 2).unwrap();
    // L[0,0] should be sqrt(4) = 2
    assert!((l[0] - 2.0).abs() < 1e-10);
}

#[test]
fn test_generalized_cholesky_non_spd_fails() {
    let m = vec![-1.0, 0.0, 0.0, 1.0];
    assert!(GeneralizedEigenvalue::cholesky(&m, 2).is_err());
}

// ── ModalAnalysis ─────────────────────────────────────────────────────────

#[test]
fn test_modal_analysis_frequency_hz() {
    let omega = 2.0 * std::f64::consts::PI; // 1 Hz
    let f = ModalAnalysis::frequency_hz(omega);
    assert!((f - 1.0).abs() < 1e-10);
}

#[test]
fn test_modal_analysis_run_counts() {
    let k = diag_mat(&[1.0, 4.0]);
    let m = identity(2);
    let cfg = EigensolverConfig {
        n_modes: 2,
        max_iter: 500,
        tolerance: 1e-10,
        ..Default::default()
    };
    let ma = ModalAnalysis::new(k, m, 2, cfg);
    let result = ma.run().unwrap();
    assert!(result.natural_frequencies.len() <= 2);
    assert_eq!(result.natural_frequencies.len(), result.mode_shapes.len());
    assert_eq!(result.natural_frequencies.len(), result.modal_masses.len());
}

#[test]
fn test_modal_analysis_nonnegative_frequencies() {
    let k = diag_mat(&[1.0, 4.0, 9.0]);
    let m = identity(3);
    let cfg = EigensolverConfig {
        n_modes: 3,
        max_iter: 500,
        tolerance: 1e-10,
        ..Default::default()
    };
    let ma = ModalAnalysis::new(k, m, 3, cfg);
    let result = ma.run().unwrap();
    for &omega in &result.natural_frequencies {
        assert!(omega >= 0.0, "frequency must be non-negative: {omega}");
    }
}

// ── FrequencyResponse ─────────────────────────────────────────────────────

#[test]
fn test_frequency_response_static_limit() {
    // At ω=0, response = K⁻¹ F  (static solution)
    let k = diag_mat(&[2.0, 8.0]);
    let m = identity(2);
    let frf = FrequencyResponse::new(k, m, 2, 0.0);
    let f = vec![2.0, 8.0];
    let u = frf.response(0.0, &f).unwrap();
    assert!((u[0] - 1.0).abs() < 1e-10);
    assert!((u[1] - 1.0).abs() < 1e-10);
}

#[test]
fn test_frequency_response_singular_at_resonance() {
    // K = diag(1), M = I → resonance at ω=1
    let k = diag_mat(&[1.0]);
    let m = identity(1);
    let frf = FrequencyResponse::new(k, m, 1, 0.0);
    let result = frf.response(1.0, &[1.0]);
    assert!(result.is_err(), "should fail at exact resonance");
}

#[test]
fn test_frequency_response_sweep_length() {
    let k = diag_mat(&[4.0, 9.0]);
    let m = identity(2);
    let frf = FrequencyResponse::new(k, m, 2, 0.0);
    let omegas: Vec<f64> = (1..=10).map(|i| i as f64 * 0.1).collect();
    let sweep = frf.sweep(&omegas, &[1.0, 1.0]);
    assert_eq!(sweep.len(), 10);
}

#[test]
fn test_frequency_response_peak_norm_positive() {
    let k = diag_mat(&[4.0, 9.0]);
    let m = identity(2);
    let frf = FrequencyResponse::new(k, m, 2, 0.0);
    let omegas: Vec<f64> = (1..=20).map(|i| i as f64 * 0.1).collect();
    let peak = frf.peak_response_norm(&omegas, &[1.0, 1.0]);
    assert!(peak > 0.0);
}

// ── SpectralNorm ──────────────────────────────────────────────────────────

#[test]
fn test_spectral_norm_identity() {
    let sn = SpectralNorm::new(identity(3), 3);
    let norm = sn.estimate(500, 1e-10);
    assert!(
        (norm - 1.0).abs() < 0.1,
        "spectral norm of I = 1, got {norm}"
    );
}

#[test]
fn test_spectral_norm_diagonal() {
    // diag(1,2,3) — spectral norm = 3
    let a = diag_mat(&[1.0, 2.0, 3.0]);
    let sn = SpectralNorm::new(a, 3);
    let norm = sn.estimate(500, 1e-10);
    assert!(
        (norm - 3.0).abs() < 0.5,
        "spectral norm should be ~3, got {norm}"
    );
}

#[test]
fn test_spectral_norm_nonnegative() {
    let (a, n) = mat2_spd();
    let sn = SpectralNorm::new(a, n);
    let norm = sn.estimate(200, 1e-8);
    assert!(norm >= 0.0);
}

// ── ModalTruncation ───────────────────────────────────────────────────────

#[test]
fn test_modal_truncation_reconstruct_zero() {
    let modes = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
    let freqs = vec![1.0, 2.0];
    let masses = vec![1.0, 1.0];
    let mt = ModalTruncation::new(modes, freqs, masses, 2);
    let u = mt.reconstruct(&[0.0, 0.0]);
    assert_eq!(u, vec![0.0, 0.0]);
}

#[test]
fn test_modal_truncation_reconstruct_first_mode() {
    let modes = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
    let freqs = vec![1.0, 2.0];
    let masses = vec![1.0, 1.0];
    let mt = ModalTruncation::new(modes, freqs, masses, 2);
    let u = mt.reconstruct(&[3.0, 0.0]);
    assert!((u[0] - 3.0).abs() < 1e-12);
    assert!(u[1].abs() < 1e-12);
}

#[test]
fn test_modal_truncation_static_response() {
    // K = diag(1,4), modes = unit, frequencies = 1,2, masses = 1,1
    let modes = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
    let freqs = vec![1.0, 2.0];
    let masses = vec![1.0, 1.0];
    let mt = ModalTruncation::new(modes, freqs, masses, 2);
    // Force [1,0] → modal force = [1,0], modal k = [1,4]
    // q = [1/1, 0/4] = [1, 0] → u = [1, 0]
    let u = mt.static_response(&[1.0, 0.0]);
    assert!((u[0] - 1.0).abs() < 1e-10);
    assert!(u[1].abs() < 1e-10);
}

#[test]
fn test_modal_truncation_mass_fraction() {
    let modes = vec![vec![1.0, 0.0]];
    let freqs = vec![1.0];
    let masses = vec![3.0];
    let mt = ModalTruncation::new(modes, freqs, masses, 2);
    let frac = mt.modal_mass_fraction(10.0);
    assert!((frac[0] - 0.3).abs() < 1e-10);
}

#[test]
fn test_modal_truncation_modal_force() {
    let modes = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
    let freqs = vec![1.0, 2.0];
    let masses = vec![1.0, 1.0];
    let mt = ModalTruncation::new(modes, freqs, masses, 2);
    let mf = mt.modal_force(&[5.0, 7.0]);
    assert!((mf[0] - 5.0).abs() < 1e-12);
    assert!((mf[1] - 7.0).abs() < 1e-12);
}

// ── Utility functions ─────────────────────────────────────────────────────

#[test]
fn test_dot_product() {
    assert!((dot(&[1.0, 2.0, 3.0], &[4.0, 5.0, 6.0]) - 32.0).abs() < 1e-12);
}

#[test]
fn test_norm2() {
    assert!((norm2(&[3.0, 4.0]) - 5.0).abs() < 1e-12);
}

#[test]
fn test_normalize() {
    let mut v = vec![3.0, 4.0];
    let n = normalize(&mut v);
    assert!((n - 5.0).abs() < 1e-12);
    assert!((norm2(&v) - 1.0).abs() < 1e-12);
}

#[test]
fn test_dense_matvec_identity() {
    let a = identity(3);
    let x = vec![1.0, 2.0, 3.0];
    let y = dense_matvec(&a, &x, 3);
    assert_eq!(y, x);
}

#[test]
fn test_lu_factor_solve() {
    // [[2,1],[4,3]] x = [1,1] → x = [1,-1]
    let mut a = vec![2.0, 1.0, 4.0, 3.0];
    let perm = lu_factor(&mut a, 2, FRF_PIVOT_TOL).unwrap();
    let x = lu_solve(&a, &perm, &[1.0, 1.0], 2);
    assert!((x[0] - 1.0).abs() < 1e-10, "x[0]={}", x[0]);
    assert!((x[1] + 1.0).abs() < 1e-10, "x[1]={}", x[1]);
}

#[test]
fn test_lu_factor_identity() {
    let mut a = identity(3);
    let res = lu_factor(&mut a, 3, FRF_PIVOT_TOL);
    assert_eq!(res, Some(vec![0, 1, 2]));
}

#[test]
fn test_tridiagonal_eigen_sorted() {
    let d = vec![2.0, 4.0, 6.0];
    let e = vec![0.5, 0.5];
    let (eigs, _) = tridiagonal_eigen(&d, &e).unwrap();
    for w in eigs.windows(2) {
        assert!(w[0] <= w[1]);
    }
}

#[test]
fn test_tridiagonal_eigen_diagonal_case() {
    // Zero off-diagonal → eigenvalues = diagonal
    let d = vec![1.0, 3.0, 2.0];
    let e = vec![0.0, 0.0];
    let (eigs, last) = tridiagonal_eigen(&d, &e).unwrap();
    assert_eq!(eigs.len(), 3);
    assert!((eigs[0] - 1.0).abs() < 1e-8);
    assert!((eigs[1] - 2.0).abs() < 1e-8);
    assert!((eigs[2] - 3.0).abs() < 1e-8);
    // The eigenvector of the value 2 is the last unit vector.
    assert!((last[1].abs() - 1.0).abs() < 1e-12);
    assert!(last[0].abs() < 1e-12 && last[2].abs() < 1e-12);
}

// ── Reference problems with known spectra ─────────────────────────────────

use std::f64::consts::PI;

/// Stiffness matrix of the fixed-fixed chain: `tridiag(-1, 2, -1)`.
fn chain_stiffness(n: usize) -> Vec<f64> {
    let mut a = vec![0.0_f64; n * n];
    for i in 0..n {
        a[i * n + i] = 2.0;
        if i + 1 < n {
            a[i * n + i + 1] = -1.0;
            a[(i + 1) * n + i] = -1.0;
        }
    }
    a
}

/// Exact eigenvalue `2 − 2cos(jπ/(n+1))` (1-based `j`) of the chain.
fn chain_eigenvalue(j: usize, n: usize) -> f64 {
    2.0 - 2.0 * (j as f64 * PI / (n as f64 + 1.0)).cos()
}

/// Linear finite-element bar of unit length with `n` free nodes:
/// `K = tridiag(-1, 2, -1)/h`, consistent `M = h/6·tridiag(1, 4, 1)`,
/// `h = 1/(n+1)`.
fn fe_bar(n: usize) -> (Vec<f64>, Vec<f64>) {
    let h = 1.0 / (n as f64 + 1.0);
    let mut k = vec![0.0_f64; n * n];
    let mut m = vec![0.0_f64; n * n];
    for i in 0..n {
        k[i * n + i] = 2.0 / h;
        m[i * n + i] = 4.0 * h / 6.0;
        if i + 1 < n {
            k[i * n + i + 1] = -1.0 / h;
            k[(i + 1) * n + i] = -1.0 / h;
            m[i * n + i + 1] = h / 6.0;
            m[(i + 1) * n + i] = h / 6.0;
        }
    }
    (k, m)
}

/// Exact eigenvalue `6/h²·(1 − cosθ)/(2 + cosθ)`, `θ = jπ/(n+1)`, of the
/// finite-element bar of [`fe_bar`].
fn fe_bar_eigenvalue(j: usize, n: usize) -> f64 {
    let h = 1.0 / (n as f64 + 1.0);
    let theta = j as f64 * PI / (n as f64 + 1.0);
    6.0 / (h * h) * (1.0 - theta.cos()) / (2.0 + theta.cos())
}

/// Deterministic stream of numbers in `[0, 1)` for reproducible test data.
fn lcg_stream(seed: u64) -> impl FnMut() -> f64 {
    let mut state = seed;
    move || {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((state >> 33) as f64) / ((1u64 << 31) as f64)
    }
}

/// Random symmetric positive-definite matrix `BᵀB + n·I`.
fn random_spd(n: usize, seed: u64) -> Vec<f64> {
    let mut next = lcg_stream(seed);
    let b: Vec<f64> = (0..n * n).map(|_| next() - 0.5).collect();
    let mut a = vec![0.0_f64; n * n];
    for i in 0..n {
        for j in 0..n {
            let mut s = 0.0;
            for p in 0..n {
                s += b[p * n + i] * b[p * n + j];
            }
            a[i * n + j] = s;
        }
        a[i * n + i] += n as f64;
    }
    a
}

/// Quadratic/bilinear form `xᵀ A y`.
fn bilinear(a: &[f64], x: &[f64], y: &[f64], n: usize) -> f64 {
    dot(x, &dense_matvec(a, y, n))
}

fn config_with_modes(n_modes: usize) -> EigensolverConfig {
    EigensolverConfig {
        n_modes,
        ..Default::default()
    }
}

// ── GeneralizedEigenvalue::solve_smallest (F1152) ─────────────────────────

#[test]
fn test_solve_smallest_matches_chain_modes_50dof() {
    let n = 50;
    let gev = GeneralizedEigenvalue::new(chain_stiffness(n), identity(n), n, config_with_modes(5));
    let (eigs, shapes) = gev.solve_smallest().unwrap();
    assert_eq!(eigs.len(), 5);
    assert_eq!(shapes.len(), 5);
    for (idx, &lambda) in eigs.iter().enumerate() {
        let exact = chain_eigenvalue(idx + 1, n);
        assert!(
            (lambda - exact).abs() < 1e-10,
            "mode {}: {lambda} vs exact {exact}",
            idx + 1
        );
    }
    // The old k-step Ritz values were [0.013, 0.405, 1.40, 2.63, 3.62].
    assert!(eigs[4] < 0.1, "fifth mode must be the fifth *smallest*");
}

#[test]
fn test_solve_smallest_shapes_satisfy_the_eigenproblem_with_consistent_mass() {
    let n = 30;
    let (k, m) = fe_bar(n);
    let gev = GeneralizedEigenvalue::new(k.clone(), m.clone(), n, config_with_modes(6));
    let (eigs, shapes) = gev.solve_smallest().unwrap();
    let k_scale = infinity_norm(&k, n);
    for (idx, (&lambda, phi)) in eigs.iter().zip(&shapes).enumerate() {
        let exact = fe_bar_eigenvalue(idx + 1, n);
        assert!(
            (lambda - exact).abs() < 1e-9 * exact,
            "mode {}: {lambda} vs exact {exact}",
            idx + 1
        );
        // K φ = λ M φ
        let kphi = dense_matvec(&k, phi, n);
        let mphi = dense_matvec(&m, phi, n);
        let residual = kphi
            .iter()
            .zip(&mphi)
            .map(|(a, b)| (a - lambda * b).abs())
            .fold(0.0_f64, f64::max);
        let phi_max = phi.iter().fold(0.0_f64, |s, v| s.max(v.abs()));
        assert!(
            residual < 1e-10 * k_scale * phi_max,
            "mode {} residual {residual}",
            idx + 1
        );
    }
    // Mass-orthonormal and stiffness-orthogonal.
    for (i, pi) in shapes.iter().enumerate() {
        for (j, pj) in shapes.iter().enumerate() {
            let mass = bilinear(&m, pi, pj, n);
            let stiff = bilinear(&k, pi, pj, n);
            let (want_m, want_k) = if i == j { (1.0, eigs[i]) } else { (0.0, 0.0) };
            assert!((mass - want_m).abs() < 1e-10, "φ{i}ᵀMφ{j} = {mass}");
            assert!(
                (stiff - want_k).abs() < 1e-8 * eigs[i].max(1.0),
                "φ{i}ᵀKφ{j} = {stiff}"
            );
        }
    }
}

#[test]
fn test_solve_smallest_keeps_multiplicity_of_degenerate_modes() {
    // K = diag(2, 2, 2, 5, 7), M = diag(1, 1, 1, 1, 1)
    let k = diag_mat(&[7.0, 2.0, 5.0, 2.0, 2.0]);
    let gev = GeneralizedEigenvalue::new(k, identity(5), 5, config_with_modes(4));
    let (eigs, shapes) = gev.solve_smallest().unwrap();
    for (got, want) in eigs.iter().zip([2.0, 2.0, 2.0, 5.0]) {
        assert!((got - want).abs() < 1e-12, "{eigs:?}");
    }
    for (i, pi) in shapes.iter().enumerate() {
        for (j, pj) in shapes.iter().enumerate() {
            let want = if i == j { 1.0 } else { 0.0 };
            assert!((dot(pi, pj) - want).abs() < 1e-12);
        }
    }
}

#[test]
fn test_solve_smallest_is_independent_of_the_unit_system() {
    let n = 20;
    let k: Vec<f64> = chain_stiffness(n).iter().map(|v| v * 1.0e12).collect();
    let m: Vec<f64> = identity(n).iter().map(|v| v * 1.0e3).collect();
    let gev = GeneralizedEigenvalue::new(k, m, n, config_with_modes(4));
    let (eigs, shapes) = gev.solve_smallest().unwrap();
    for (idx, &lambda) in eigs.iter().enumerate() {
        let exact = 1.0e9 * chain_eigenvalue(idx + 1, n);
        assert!(
            (lambda - exact).abs() < 1e-9 * exact,
            "mode {}: {lambda} vs {exact}",
            idx + 1
        );
    }
    // Unit modal mass: φᵀ(1e3·I)φ = 1.
    let mass: f64 = shapes[0].iter().map(|v| v * v * 1.0e3).sum();
    assert!((mass - 1.0).abs() < 1e-12, "modal mass {mass}");
}

#[test]
fn test_solve_smallest_sign_convention_and_trivial_requests() {
    let n = 12;
    let gev = GeneralizedEigenvalue::new(chain_stiffness(n), identity(n), n, config_with_modes(12));
    let (eigs, shapes) = gev.solve_smallest().unwrap();
    assert_eq!(eigs.len(), n);
    for phi in &shapes {
        let dominant = phi
            .iter()
            .fold(0.0_f64, |b, &v| if v.abs() > b.abs() { v } else { b });
        assert!(dominant > 0.0, "dominant component must be positive");
    }
    // More modes than DOFs are clamped, zero modes give nothing.
    let gev = GeneralizedEigenvalue::new(chain_stiffness(3), identity(3), 3, config_with_modes(10));
    assert_eq!(gev.solve_smallest().unwrap().0.len(), 3);
    let gev = GeneralizedEigenvalue::new(chain_stiffness(3), identity(3), 3, config_with_modes(0));
    let (e, s) = gev.solve_smallest().unwrap();
    assert!(e.is_empty() && s.is_empty());
}

#[test]
fn test_solve_smallest_rejects_invalid_matrices() {
    let cfg = config_with_modes(2);
    // Asymmetric stiffness.
    let k = vec![2.0, -1.0, -0.5, 2.0];
    let err = GeneralizedEigenvalue::new(k, identity(2), 2, cfg.clone())
        .solve_smallest()
        .unwrap_err();
    assert!(err.contains("not symmetric"), "{err}");
    // Asymmetric mass.
    let m = vec![2.0, 0.5, 0.0, 2.0];
    let err = GeneralizedEigenvalue::new(identity(2), m, 2, cfg.clone())
        .solve_smallest()
        .unwrap_err();
    assert!(err.contains("not symmetric"), "{err}");
    // Non-finite entries.
    let k = vec![f64::NAN, 0.0, 0.0, 1.0];
    let err = GeneralizedEigenvalue::new(k, identity(2), 2, cfg.clone())
        .solve_smallest()
        .unwrap_err();
    assert!(err.contains("non-finite"), "{err}");
    // Mass matrix that is not positive definite.
    let m = vec![1.0, 2.0, 2.0, 1.0];
    let err = GeneralizedEigenvalue::new(identity(2), m, 2, cfg)
        .solve_smallest()
        .unwrap_err();
    assert!(err.contains("positive definite"), "{err}");
}

// ── ModalAnalysis (F1152) ─────────────────────────────────────────────────

#[test]
fn test_modal_analysis_reports_the_lowest_chain_frequencies() {
    let n = 50;
    let ma = ModalAnalysis::new(chain_stiffness(n), identity(n), n, config_with_modes(5));
    let result = ma.run().unwrap();
    assert_eq!(result.natural_frequencies.len(), 5);
    for (idx, &omega) in result.natural_frequencies.iter().enumerate() {
        let exact = chain_eigenvalue(idx + 1, n).sqrt();
        assert!(
            (omega - exact).abs() < 1e-10,
            "ω{}: {omega} vs {exact}",
            idx + 1
        );
    }
    for (idx, (&mm, &km)) in result
        .modal_masses
        .iter()
        .zip(&result.modal_stiffnesses)
        .enumerate()
    {
        assert!((mm - 1.0).abs() < 1e-10, "modal mass {idx}: {mm}");
        let lambda = chain_eigenvalue(idx + 1, n);
        assert!(
            (km - lambda).abs() < 1e-9 * lambda.max(1e-3),
            "modal k {idx}: {km}"
        );
    }
}

#[test]
fn test_modal_analysis_participation_factors_are_complete() {
    // With all modes retained, Σ Γᵢ² mᵢ = rᵀ M r  (completeness of the modes).
    let n = 8;
    let (k, m) = fe_bar(n);
    let ma = ModalAnalysis::new(k, m.clone(), n, config_with_modes(n));
    let result = ma.run().unwrap();
    let total: f64 = result
        .participation_factors
        .iter()
        .zip(&result.modal_masses)
        .map(|(g, mm)| g * g * mm)
        .sum();
    let ones = vec![1.0_f64; n];
    let expected = bilinear(&m, &ones, &ones, n);
    assert!(
        (total - expected).abs() < 1e-10 * expected,
        "Σ Γ² m = {total}, rᵀMr = {expected}"
    );
}

// ── LanczosMethod (F1152) ─────────────────────────────────────────────────

#[test]
fn test_lanczos_recovers_the_smallest_chain_eigenvalues() {
    let n = 50;
    let lm = LanczosMethod::new(chain_stiffness(n), n, config_with_modes(5));
    let eigs = lm.eigenvalues();
    assert_eq!(eigs.len(), 5);
    for (idx, &lambda) in eigs.iter().enumerate() {
        let exact = chain_eigenvalue(idx + 1, n);
        assert!(
            (lambda - exact).abs() < 1e-9,
            "λ{}: {lambda} vs exact {exact}",
            idx + 1
        );
    }
}

#[test]
fn test_lanczos_matches_dense_solver_on_a_random_spd_matrix() {
    let n = 40;
    let a = random_spd(n, 0x5EED);
    let k = 6;
    let lm = LanczosMethod::new(a.clone(), n, config_with_modes(k));
    let eigs = lm.eigenvalues();
    let dense = GeneralizedEigenvalue::new(a, identity(n), n, config_with_modes(k));
    let (reference, _) = dense.solve_smallest().unwrap();
    assert_eq!(eigs.len(), k);
    for (got, want) in eigs.iter().zip(&reference) {
        assert!(
            (got - want).abs() < 1e-9 * want.abs().max(1.0),
            "{got} vs {want}"
        );
    }
}

#[test]
fn test_lanczos_finds_isolated_smallest_eigenvalues_of_a_large_matrix() {
    // Three small eigenvalues far below a tight cluster of 147 others.
    let n = 150;
    let mut diag: Vec<f64> = (0..n).map(|i| 100.0 + 0.01 * i as f64).collect();
    diag[17] = 1.0;
    diag[90] = 2.5;
    diag[3] = 4.0;
    let lm = LanczosMethod::new(diag_mat(&diag), n, config_with_modes(3));
    let eigs = lm.eigenvalues();
    for (got, want) in eigs.iter().zip([1.0, 2.5, 4.0]) {
        assert!((got - want).abs() < 1e-9, "{eigs:?}");
    }
}

#[test]
fn test_lanczos_restarts_after_invariant_subspace_to_keep_multiplicity() {
    // The identity has a one-dimensional Krylov space; the recurrence must
    // restart to deliver three values.
    let eigs = LanczosMethod::new(identity(4), 4, config_with_modes(3)).eigenvalues();
    assert_eq!(eigs.len(), 3);
    for e in &eigs {
        assert!((e - 1.0).abs() < 1e-12, "{eigs:?}");
    }
    // diag(1,1,1,2,2,3): three distinct values, but k = 4 needs the repeats.
    let a = diag_mat(&[3.0, 2.0, 1.0, 2.0, 1.0, 1.0]);
    let eigs = LanczosMethod::new(a, 6, config_with_modes(4)).eigenvalues();
    for (got, want) in eigs.iter().zip([1.0, 1.0, 1.0, 2.0]) {
        assert!((got - want).abs() < 1e-12, "{eigs:?}");
    }
}

#[test]
fn test_lanczos_handles_indefinite_and_degenerate_inputs() {
    let a = diag_mat(&[5.0, -3.0, 2.0, -1.0, 7.0]);
    let eigs = LanczosMethod::new(a, 5, config_with_modes(2)).eigenvalues();
    assert!(
        (eigs[0] + 3.0).abs() < 1e-12 && (eigs[1] + 1.0).abs() < 1e-12,
        "{eigs:?}"
    );
    // More modes requested than the matrix has.
    let eigs = LanczosMethod::new(chain_stiffness(3), 3, config_with_modes(9)).eigenvalues();
    assert_eq!(eigs.len(), 3);
    // Zero modes, the zero matrix and non-finite input.
    assert!(
        LanczosMethod::new(identity(3), 3, config_with_modes(0))
            .eigenvalues()
            .is_empty()
    );
    let zero = vec![0.0_f64; 9];
    assert_eq!(
        LanczosMethod::new(zero, 3, config_with_modes(2)).eigenvalues(),
        vec![0.0, 0.0]
    );
    let mut bad = identity(3);
    bad[4] = f64::NAN;
    let eigs = LanczosMethod::new(bad, 3, config_with_modes(2)).eigenvalues();
    assert!(eigs.len() == 2 && eigs.iter().all(|e| e.is_nan()));
}

// ── Tridiagonal kernel ────────────────────────────────────────────────────

#[test]
fn test_tridiagonal_eigen_chain_values_and_last_components() {
    let n = 50;
    let diag = vec![2.0_f64; n];
    let off = vec![-1.0_f64; n - 1];
    let (values, last) = tridiagonal_eigen(&diag, &off).unwrap();
    for (j, (&lambda, &s)) in values.iter().zip(&last).enumerate() {
        let idx = j + 1;
        assert!((lambda - chain_eigenvalue(idx, n)).abs() < 1e-12, "λ{idx}");
        // Last component of the normalised eigenvector sqrt(2/(n+1))·sin(jπ/(n+1)).
        let exact = (2.0 / (n as f64 + 1.0)).sqrt() * (idx as f64 * PI / (n as f64 + 1.0)).sin();
        assert!((s.abs() - exact).abs() < 1e-12, "last component {idx}: {s}");
    }
}

#[test]
fn test_ritz_pairs_jacobi_fallback_agrees_with_ql() {
    let mut next = lcg_stream(77);
    let m = 12;
    let alpha: Vec<f64> = (0..m).map(|_| 4.0 * next() - 1.0).collect();
    let beta: Vec<f64> = (0..m - 1).map(|_| 0.2 + next()).collect();
    let (ql_values, ql_last) = tridiagonal_eigen(&alpha, &beta).unwrap();
    let (jac_values, jac_last) = ritz_pairs_jacobi(&alpha, &beta);
    for i in 0..m {
        assert!((ql_values[i] - jac_values[i]).abs() < 1e-11, "value {i}");
        assert!(
            (ql_last[i].abs() - jac_last[i].abs()).abs() < 1e-9,
            "last {i}"
        );
    }
}

#[test]
fn test_tridiagonal_eigen_splits_blocks_and_handles_trivial_sizes() {
    assert_eq!(tridiagonal_eigen(&[], &[]), Some((vec![], vec![])));
    let (v, last) = tridiagonal_eigen(&[3.5], &[]).unwrap();
    assert_eq!((v, last), (vec![3.5], vec![1.0]));
    // A zero coupling decouples the matrix into [[2,1],[1,2]] ⊕ [5].
    let (values, _) = tridiagonal_eigen(&[2.0, 2.0, 5.0], &[1.0, 0.0]).unwrap();
    for (got, want) in values.iter().zip([1.0, 3.0, 5.0]) {
        assert!((got - want).abs() < 1e-13, "{values:?}");
    }
}

// ── Pivoted LU and FrequencyResponse (F1155) ──────────────────────────────

#[test]
fn test_lu_factor_pivots_around_a_zero_diagonal() {
    // [[0,1],[1,0]] x = [3,5]  →  x = [5,3]; the unpivoted LU failed here.
    let mut a = vec![0.0, 1.0, 1.0, 0.0];
    let perm = lu_factor(&mut a, 2, FRF_PIVOT_TOL).expect("nonsingular");
    assert_eq!(perm, vec![1, 0]);
    let x = lu_solve(&a, &perm, &[3.0, 5.0], 2);
    assert_eq!(x, vec![5.0, 3.0]);
}

#[test]
fn test_lu_solve_residual_on_a_matrix_that_needs_many_swaps() {
    let n = 5;
    let mut next = lcg_stream(4242);
    let a: Vec<f64> = (0..n * n).map(|_| next() - 0.5).collect();
    let b: Vec<f64> = (0..n).map(|i| (i as f64 + 1.0) * 0.7).collect();
    let mut lu = a.clone();
    let perm = lu_factor(&mut lu, n, FRF_PIVOT_TOL).expect("random matrix is nonsingular");
    let x = lu_solve(&lu, &perm, &b, n);
    let ax = dense_matvec(&a, &x, n);
    for (got, want) in ax.iter().zip(&b) {
        assert!((got - want).abs() < 1e-12, "A x = {ax:?}");
    }
}

#[test]
fn test_lu_factor_rejects_singular_and_non_finite_matrices() {
    let mut rank_one = vec![1.0, 2.0, 2.0, 4.0];
    assert_eq!(lu_factor(&mut rank_one, 2, FRF_PIVOT_TOL), None);
    // An all-zero matrix has a zero threshold; the zero pivot must still fail.
    let mut zero = vec![0.0_f64; 4];
    assert_eq!(lu_factor(&mut zero, 2, FRF_PIVOT_TOL), None);
    assert_eq!(lu_factor(&mut zero, 2, 0.0), None);
    let mut nan = vec![1.0, f64::NAN, 0.0, 1.0];
    assert_eq!(lu_factor(&mut nan, 2, 0.0), None);
    // The relative threshold scales with the matrix: tiny but regular is fine.
    let mut tiny = vec![1.0e-200, 0.0, 0.0, 2.0e-200];
    assert!(lu_factor(&mut tiny, 2, FRF_PIVOT_TOL).is_some());
}

#[test]
fn test_frequency_response_applies_the_damping_ratio() {
    // k = m = 1, ζ = 0.05, ω = 1: u = F/(i·2ζ) → |u| = 10, purely imaginary.
    let frf = FrequencyResponse::new(vec![1.0], vec![1.0], 1, 0.05);
    let u = frf.response(1.0, &[1.0]).unwrap();
    assert!((u[0] - 10.0).abs() < 1e-12, "|u| = {}", u[0]);
    let (re, im) = frf.response_complex(1.0, &[1.0]).unwrap();
    assert!(re[0].abs() < 1e-12, "re = {}", re[0]);
    assert!((im[0] + 10.0).abs() < 1e-12, "im = {}", im[0]);
    assert!((frf.loss_factor() - 0.1).abs() < 1e-15);
    // Undamped: the same frequency is a true resonance.
    let undamped = FrequencyResponse::new(vec![1.0], vec![1.0], 1, 0.0);
    assert!(undamped.response(1.0, &[1.0]).is_err());
}

#[test]
fn test_frequency_response_phase_and_off_resonance_values() {
    // u = F / (k(1 + iη) − ω²) with k = m = 1, ζ = 0.05, ω = 0.5.
    let frf = FrequencyResponse::new(vec![1.0], vec![1.0], 1, 0.05);
    let (re, im) = frf.response_complex(0.5, &[2.0]).unwrap();
    let (dr, di) = (0.75_f64, 0.1_f64);
    let den = dr * dr + di * di;
    assert!((re[0] - 2.0 * dr / den).abs() < 1e-12, "re = {}", re[0]);
    assert!((im[0] + 2.0 * di / den).abs() < 1e-12, "im = {}", im[0]);
    assert!(im[0] < 0.0, "response lags the load");
    let mag = frf.response(0.5, &[2.0]).unwrap();
    assert!((mag[0] - 2.0 / den.sqrt()).abs() < 1e-12);
}

#[test]
fn test_frequency_response_matches_modal_superposition_with_damping() {
    let n = 6;
    let (k, m) = fe_bar(n);
    let zeta = 0.03;
    let omega = 0.9 * fe_bar_eigenvalue(2, n).sqrt();
    let force = vec![1.0, 0.0, -2.0, 0.5, 0.0, 1.5];
    let frf = FrequencyResponse::new(k.clone(), m.clone(), n, zeta);
    let (re, im) = frf.response_complex(omega, &force).unwrap();

    let gev = GeneralizedEigenvalue::new(k, m, n, config_with_modes(n));
    let (lambdas, shapes) = gev.solve_smallest().unwrap();
    let eta = 2.0 * zeta;
    let mut re_ref = vec![0.0_f64; n];
    let mut im_ref = vec![0.0_f64; n];
    for (lambda, phi) in lambdas.iter().zip(&shapes) {
        let load = dot(phi, &force);
        let (a, b) = (lambda - omega * omega, eta * lambda);
        let den = a * a + b * b;
        let (q_re, q_im) = (load * a / den, -load * b / den);
        for i in 0..n {
            re_ref[i] += q_re * phi[i];
            im_ref[i] += q_im * phi[i];
        }
    }
    let scale = re_ref
        .iter()
        .chain(&im_ref)
        .fold(0.0_f64, |s, v| s.max(v.abs()));
    for i in 0..n {
        assert!((re[i] - re_ref[i]).abs() < 1e-9 * scale, "re[{i}]");
        assert!((im[i] - im_ref[i]).abs() < 1e-9 * scale, "im[{i}]");
    }
}

#[test]
fn test_frequency_response_needs_pivoting_when_a_leading_pivot_vanishes() {
    // K = [[2,-1],[-1,2]], M = I, ω² = 2: K − ω²M = [[0,-1],[-1,0]] is
    // nonsingular (the eigenvalues of K are 1 and 3) but has a zero pivot.
    let k = vec![2.0, -1.0, -1.0, 2.0];
    let frf = FrequencyResponse::new(k, identity(2), 2, 0.0);
    let omega = 2.0_f64.sqrt();
    let (re, im) = frf.response_complex(omega, &[1.0, 0.0]).unwrap();
    assert!(re[0].abs() < 1e-12 && (re[1] + 1.0).abs() < 1e-12, "{re:?}");
    assert!(im.iter().all(|v| *v == 0.0));
    // Above the first resonance the magnitude is still the modulus.
    let u = frf.response(omega, &[1.0, 0.0]).unwrap();
    assert!(u[0].abs() < 1e-12 && (u[1] - 1.0).abs() < 1e-12);
}

#[test]
fn test_frequency_response_damped_peak_is_bounded_by_the_loss_factor() {
    let zeta = 0.05;
    let frf = FrequencyResponse::new(vec![1.0], vec![1.0], 1, zeta);
    let omegas: Vec<f64> = (0..=1000).map(|i| 0.5 + 0.001 * i as f64).collect();
    let peak = frf.peak_response_norm(&omegas, &[1.0]);
    assert!((peak - 1.0 / (2.0 * zeta)).abs() < 1e-2, "peak = {peak}");
    // The undamped system has no finite peak at the sampled resonance.
    let undamped = FrequencyResponse::new(vec![1.0], vec![1.0], 1, 0.0);
    let exact_resonance = undamped.sweep(&[1.0], &[1.0]);
    assert!(exact_resonance[0].1.is_err());
}

#[test]
fn test_frequency_response_rejects_invalid_arguments() {
    let frf = FrequencyResponse::new(vec![1.0], vec![1.0], 1, -0.1);
    assert!(frf.response(0.5, &[1.0]).unwrap_err().contains("damping"));
    let frf = FrequencyResponse::new(vec![1.0], vec![1.0], 1, f64::NAN);
    assert!(frf.response(0.5, &[1.0]).is_err());
    let frf = FrequencyResponse::new(vec![1.0], vec![1.0], 1, 0.02);
    assert!(frf.response(f64::NAN, &[1.0]).is_err());
    assert!(frf.response(0.5, &[1.0, 2.0]).is_err());
    let frf = FrequencyResponse::new(vec![1.0, 0.0], vec![1.0], 1, 0.02);
    assert!(frf.response(0.5, &[1.0]).is_err());
}

#[test]
fn test_inverse_iteration_survives_a_zero_leading_pivot() {
    // A = [[2,1],[1,3]], σ = 2: A − σI = [[0,1],[1,1]] has a zero pivot.
    let a = vec![2.0, 1.0, 1.0, 3.0];
    let cfg = EigensolverConfig {
        max_iter: 500,
        tolerance: 1e-14,
        shift: 2.0,
        n_modes: 1,
    };
    let (lam, v) = InverseIteration::new(a, 2, cfg).solve().unwrap();
    let exact = (5.0 - 5.0_f64.sqrt()) / 2.0;
    assert!((lam - exact).abs() < 1e-10, "λ = {lam}, exact {exact}");
    assert!((norm2(&v) - 1.0).abs() < 1e-12);
}
