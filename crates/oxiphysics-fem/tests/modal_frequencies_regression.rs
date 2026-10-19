// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for [`ModalAnalysis`].
//!
//! The previous Rayleigh-quotient iteration shifted by the quotient of a sine
//! start vector, so it converged to whichever mode happened to be nearest: for a
//! 10-DOF spring chain asked for three modes it returned the unsorted
//! `[0.831, 1.081, 0.563]` instead of `[0.285, 0.563, 0.831]` (fundamental
//! missing).  `compute_mode_shapes` returned its *start vector* whenever
//! `K − ω²·M` was singular, e.g. exactly at a natural frequency.  The tests
//! compare against closed-form spectra and the defining eigen-relations.

use oxiphysics_fem::analysis::{FemError, ModalAnalysis};
use std::f64::consts::PI;

/// `ks · tridiag(−1, 2, −1)`: a fixed-fixed chain of `n` unit springs.
fn spring_chain_stiffness(n: usize, ks: f64) -> Vec<f64> {
    let mut k = vec![0.0; n * n];
    for i in 0..n {
        k[i * n + i] = 2.0 * ks;
        if i + 1 < n {
            k[i * n + i + 1] = -ks;
            k[(i + 1) * n + i] = -ks;
        }
    }
    k
}

fn scaled_identity(n: usize, scale: f64) -> Vec<f64> {
    let mut m = vec![0.0; n * n];
    for i in 0..n {
        m[i * n + i] = scale;
    }
    m
}

fn matvec(a: &[f64], x: &[f64]) -> Vec<f64> {
    let n = x.len();
    (0..n)
        .map(|i| (0..n).map(|j| a[i * n + j] * x[j]).sum())
        .collect()
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// The component that fixes the sign of a mode shape: the first one whose
/// magnitude is within `1e-9` (relative) of the largest.
fn sign_reference(shape: &[f64]) -> f64 {
    let largest = shape.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    shape
        .iter()
        .copied()
        .find(|v| v.abs() >= (1.0 - 1.0e-9) * largest)
        .unwrap_or(0.0)
}

/// Exact `ω_j = 2·sqrt(ks/m)·sin(jπ / (2(n+1)))` of the fixed-fixed chain.
fn chain_frequency(j: usize, n: usize, ks: f64, mass: f64) -> f64 {
    2.0 * (ks / mass).sqrt() * (j as f64 * PI / (2.0 * (n as f64 + 1.0))).sin()
}

#[test]
fn lowest_chain_frequencies_are_exact_and_ascending_without_sorting() {
    let (n, ks, mass) = (10, 3.0e6, 2.5);
    let k = spring_chain_stiffness(n, ks);
    let m = scaled_identity(n, mass);
    let freqs = ModalAnalysis::compute_natural_frequencies(&k, &m, n, 3, 500)
        .expect("SPD mass matrix, symmetric stiffness");
    assert_eq!(freqs.len(), 3);
    for (index, &omega) in freqs.iter().enumerate() {
        let expected = chain_frequency(index + 1, n, ks, mass);
        assert!(
            (omega - expected).abs() <= 1e-10 * expected,
            "mode {}: ω = {omega}, exact {expected}",
            index + 1
        );
    }
    assert!(
        freqs.windows(2).all(|pair| pair[0] < pair[1]),
        "frequencies must be returned in ascending order: {freqs:?}"
    );
}

#[test]
fn every_mode_of_the_chain_is_found_exactly_once() {
    let (n, ks, mass) = (10, 1.0, 1.0);
    let k = spring_chain_stiffness(n, ks);
    let m = scaled_identity(n, mass);
    let freqs = ModalAnalysis::compute_natural_frequencies(&k, &m, n, n, 500).expect("modal");
    assert_eq!(freqs.len(), n);
    for (index, &omega) in freqs.iter().enumerate() {
        let expected = chain_frequency(index + 1, n, ks, mass);
        assert!((omega - expected).abs() <= 1e-10 * expected, "mode {index}");
    }
}

#[test]
fn requesting_more_modes_than_dofs_returns_all_dofs() {
    let k = spring_chain_stiffness(4, 1.0);
    let m = scaled_identity(4, 1.0);
    let freqs = ModalAnalysis::compute_natural_frequencies(&k, &m, 4, 9, 200).expect("modal");
    assert_eq!(freqs.len(), 4);
    let none = ModalAnalysis::compute_natural_frequencies(&k, &m, 4, 0, 200).expect("modal");
    assert!(none.is_empty());
}

#[test]
fn physical_si_units_converge_as_fast_as_unit_matrices() {
    // Stiffness ~ 1e9 N/m, mass ~ 1e-3 kg: an absolute Jacobi threshold would
    // never be reached and the old absolute pivot test mis-judged the scale.
    let (n, ks, mass) = (12, 4.0e9, 3.0e-3);
    let k = spring_chain_stiffness(n, ks);
    let m = scaled_identity(n, mass);
    let freqs = ModalAnalysis::compute_natural_frequencies(&k, &m, n, 4, 50).expect("modal");
    for (index, &omega) in freqs.iter().enumerate() {
        let expected = chain_frequency(index + 1, n, ks, mass);
        assert!((omega - expected).abs() <= 1e-9 * expected, "mode {index}");
    }
}

#[test]
fn non_diagonal_pencil_matches_the_hand_computed_spectrum() {
    // det(K − λM) = (3 − 2λ)(1 − λ) − 1 = 2λ² − 5λ + 2 ⇒ λ ∈ {1/2, 2}.
    let k = [3.0, -1.0, -1.0, 1.0];
    let m = [2.0, 0.0, 0.0, 1.0];
    let freqs = ModalAnalysis::compute_natural_frequencies(&k, &m, 2, 2, 100).expect("modal");
    assert!((freqs[0] - 0.5_f64.sqrt()).abs() < 1e-14, "{freqs:?}");
    assert!((freqs[1] - 2.0_f64.sqrt()).abs() < 1e-14, "{freqs:?}");
}

#[test]
fn consistent_mass_free_bar_has_a_rigid_mode_and_the_textbook_frequency() {
    // One bar element, consistent mass: ω² = 12·k/m (m = ρAL), plus a rigid mode.
    let (k_axial, total_mass) = (7.0e4, 3.0);
    let k = [k_axial, -k_axial, -k_axial, k_axial];
    let m = [
        total_mass / 3.0,
        total_mass / 6.0,
        total_mass / 6.0,
        total_mass / 3.0,
    ];
    let freqs = ModalAnalysis::compute_natural_frequencies(&k, &m, 2, 2, 100).expect("modal");
    let expected = (12.0 * k_axial / total_mass).sqrt();
    assert!((freqs[1] - expected).abs() <= 1e-12 * expected, "{freqs:?}");
    // The rigid-body mode is zero up to the square root of round-off.
    assert!(
        freqs[0] <= 1e-6 * freqs[1],
        "rigid-body mode ω = {}",
        freqs[0]
    );
}

#[test]
fn mode_shapes_are_mass_orthonormal_and_solve_the_eigenproblem() {
    let n = 6;
    let k = spring_chain_stiffness(n, 5.0e3);
    // A non-diagonal, SPD mass matrix.
    let mut m = scaled_identity(n, 2.0);
    for i in 0..n - 1 {
        m[i * n + i + 1] = 0.3;
        m[(i + 1) * n + i] = 0.3;
    }
    let modes = ModalAnalysis::compute_modes(&k, &m, n, n, 500).expect("modal");
    assert_eq!(modes.len(), n);
    for (i, (omega_i, phi_i)) in modes.iter().enumerate() {
        let kphi = matvec(&k, phi_i);
        let mphi = matvec(&m, phi_i);
        let residual: f64 = kphi
            .iter()
            .zip(&mphi)
            .map(|(a, b)| (a - omega_i * omega_i * b).powi(2))
            .sum::<f64>()
            .sqrt();
        let scale = dot(&kphi, &kphi).sqrt();
        assert!(residual <= 1e-10 * scale, "mode {i}: residual {residual:e}");
        for (j, (_, phi_j)) in modes.iter().enumerate() {
            let expected = if i == j { 1.0 } else { 0.0 };
            let inner = dot(phi_j, &mphi);
            assert!(
                (inner - expected).abs() < 1e-12,
                "φ{j}ᵀ·M·φ{i} = {inner}, expected {expected}"
            );
        }
        // Deterministic sign: the first component of (nearly) largest magnitude
        // is positive.
        assert!(
            sign_reference(phi_i) > 0.0,
            "mode {i} sign is not normalised: {phi_i:?}"
        );
    }
    assert!(modes.windows(2).all(|pair| pair[0].0 <= pair[1].0));
}

#[test]
fn repeated_eigenvalues_keep_orthonormal_shapes() {
    let n = 3;
    let identity = scaled_identity(n, 1.0);
    let modes = ModalAnalysis::compute_modes(&identity, &identity, n, n, 100).expect("modal");
    for (omega, _) in &modes {
        assert!((omega - 1.0).abs() < 1e-14);
    }
    for i in 0..n {
        for j in 0..n {
            let expected = if i == j { 1.0 } else { 0.0 };
            assert!((dot(&modes[i].1, &modes[j].1) - expected).abs() < 1e-14);
        }
    }
}

#[test]
fn mass_matrix_must_be_symmetric_positive_definite() {
    let k = spring_chain_stiffness(2, 1.0);
    let indefinite = [1.0, 2.0, 2.0, 1.0];
    assert!(matches!(
        ModalAnalysis::compute_natural_frequencies(&k, &indefinite, 2, 2, 100),
        Err(FemError::InvalidInput(_))
    ));
    // A massless DOF.
    let massless = [1.0, 0.0, 0.0, 0.0];
    assert!(matches!(
        ModalAnalysis::compute_natural_frequencies(&k, &massless, 2, 2, 100),
        Err(FemError::InvalidInput(_))
    ));
    // A non-symmetric mass matrix.
    let skew = [1.0, 0.5, 0.0, 1.0];
    assert!(matches!(
        ModalAnalysis::compute_natural_frequencies(&k, &skew, 2, 2, 100),
        Err(FemError::InvalidInput(_))
    ));
}

#[test]
fn malformed_stiffness_or_sizes_are_rejected() {
    let m = scaled_identity(2, 1.0);
    let skew_k = [2.0, -1.0, 0.0, 2.0];
    assert!(matches!(
        ModalAnalysis::compute_natural_frequencies(&skew_k, &m, 2, 1, 100),
        Err(FemError::InvalidInput(_))
    ));
    assert!(matches!(
        ModalAnalysis::compute_natural_frequencies(&[1.0; 3], &m, 2, 1, 100),
        Err(FemError::InvalidInput(_))
    ));
    let nan_k = [f64::NAN, 0.0, 0.0, 1.0];
    assert!(matches!(
        ModalAnalysis::compute_natural_frequencies(&nan_k, &m, 2, 1, 100),
        Err(FemError::InvalidInput(_))
    ));
}

#[test]
fn an_unstable_structure_has_no_real_frequency() {
    // K has a genuinely negative eigenvalue: clamping it to ω = 0 would hide the
    // instability, so it is an error.
    let k = [-1.0, 0.0, 0.0, 4.0];
    let m = scaled_identity(2, 1.0);
    assert!(matches!(
        ModalAnalysis::compute_natural_frequencies(&k, &m, 2, 2, 100),
        Err(FemError::InvalidInput(_))
    ));
}

#[test]
fn too_few_sweeps_report_non_convergence_instead_of_returning_a_partial_answer() {
    let n = 10;
    let k = spring_chain_stiffness(n, 1.0);
    let m = scaled_identity(n, 1.0);
    assert!(matches!(
        ModalAnalysis::compute_natural_frequencies(&k, &m, n, 3, 1),
        Err(FemError::SolverNotConverged)
    ));
}

#[test]
fn instance_wrappers_use_the_instance_settings() {
    let n = 5;
    let k = spring_chain_stiffness(n, 2.0);
    let m = scaled_identity(n, 1.0);
    let analysis = ModalAnalysis::default();
    let via_instance = analysis.frequencies(&k, &m, n, 3).expect("modal");
    let via_static =
        ModalAnalysis::compute_natural_frequencies(&k, &m, n, 3, analysis.max_iter).expect("modal");
    assert_eq!(via_instance, via_static);

    let omega = via_instance[1];
    let shape = analysis.mode_shape(&k, &m, n, omega).expect("mode shape");
    let direct = ModalAnalysis::compute_mode_shapes(&k, &m, n, omega, analysis.max_iter)
        .expect("mode shape");
    for (a, b) in shape.iter().zip(&direct) {
        assert!((a - b).abs() < 1e-9);
    }
}

/// `K = diag(1, 9)`, `M = I`, `ω = 1` sits exactly on the first eigenvalue, so
/// `K − ω²·M = diag(0, 8)` is singular.  The old code returned its start vector
/// `[0.495, 0.869]`; the shape is `[1, 0]`.
#[test]
fn mode_shape_at_an_exact_natural_frequency_is_the_eigenvector_not_the_start_vector() {
    let k = [1.0, 0.0, 0.0, 9.0];
    let m = [1.0, 0.0, 0.0, 1.0];
    let phi = ModalAnalysis::compute_mode_shapes(&k, &m, 2, 1.0, 100).expect("mode shape");
    assert!((phi[0] - 1.0).abs() < 1e-8, "{phi:?}");
    assert!(phi[1].abs() < 1e-8, "{phi:?}");

    let phi2 = ModalAnalysis::compute_mode_shapes(&k, &m, 2, 3.0, 100).expect("mode shape");
    assert!(phi2[0].abs() < 1e-8, "{phi2:?}");
    assert!((phi2[1] - 1.0).abs() < 1e-8, "{phi2:?}");
}

#[test]
fn mode_shape_converges_to_the_nearest_mode_for_an_inexact_frequency() {
    // ω² = 4 is closer to λ = 1 (distance 3) than to λ = 9 (distance 5).
    let k = [1.0, 0.0, 0.0, 9.0];
    let m = [1.0, 0.0, 0.0, 1.0];
    let phi = ModalAnalysis::compute_mode_shapes(&k, &m, 2, 2.0, 500).expect("mode shape");
    assert!((phi[0] - 1.0).abs() < 1e-8, "{phi:?}");
    assert!(phi[1].abs() < 1e-8, "{phi:?}");
}

#[test]
fn mode_shape_of_a_non_diagonal_pencil_is_mass_normalised() {
    // λ = 2 has the eigenvector ∝ [1, −1]; φᵀMφ = 2a² + a² = 1 ⇒ a = 1/√3.
    let k = [3.0, -1.0, -1.0, 1.0];
    let m = [2.0, 0.0, 0.0, 1.0];
    let phi = ModalAnalysis::compute_mode_shapes(&k, &m, 2, 2.0_f64.sqrt(), 100).expect("shape");
    let a = 1.0 / 3.0_f64.sqrt();
    assert!((phi[0] - a).abs() < 1e-8, "{phi:?}");
    assert!((phi[1] + a).abs() < 1e-8, "{phi:?}");
    let mphi = matvec(&m, &phi);
    assert!((dot(&phi, &mphi) - 1.0).abs() < 1e-12);
}

#[test]
fn mode_shape_agrees_with_the_direct_eigen_solution_for_a_chain() {
    let n = 8;
    let k = spring_chain_stiffness(n, 100.0);
    let m = scaled_identity(n, 0.5);
    let modes = ModalAnalysis::compute_modes(&k, &m, n, n, 500).expect("modal");
    for (omega, direct) in modes.iter().take(4) {
        let shape = ModalAnalysis::compute_mode_shapes(&k, &m, n, *omega, 500).expect("shape");
        for (a, b) in shape.iter().zip(direct) {
            assert!((a - b).abs() < 1e-6, "{shape:?} vs {direct:?}");
        }
    }
}

#[test]
fn mode_shape_reports_unsolvable_or_ill_posed_problems() {
    let k = [1.0, 0.0, 0.0, 9.0];
    let m = [1.0, 0.0, 0.0, 1.0];
    // ω² = 5 is equidistant from both eigenvalues: inverse iteration oscillates.
    assert!(matches!(
        ModalAnalysis::compute_mode_shapes(&k, &m, 2, 5.0_f64.sqrt(), 200),
        Err(FemError::SolverNotConverged)
    ));
    // A singular pencil: no shift makes K − ω²·M regular.
    let zero = [0.0; 4];
    assert!(matches!(
        ModalAnalysis::compute_mode_shapes(&zero, &zero, 2, 1.0, 50),
        Err(FemError::SingularSystem)
    ));
    // Invalid arguments.
    assert!(matches!(
        ModalAnalysis::compute_mode_shapes(&k, &m, 2, 1.0, 0),
        Err(FemError::InvalidInput(_))
    ));
    assert!(matches!(
        ModalAnalysis::compute_mode_shapes(&k, &m, 2, f64::NAN, 10),
        Err(FemError::InvalidInput(_))
    ));
    assert!(matches!(
        ModalAnalysis::compute_mode_shapes(&k, &[1.0; 3], 2, 1.0, 10),
        Err(FemError::InvalidInput(_))
    ));
    assert!(
        ModalAnalysis::compute_mode_shapes(&[], &[], 0, 1.0, 10)
            .expect("empty system")
            .is_empty()
    );
}
