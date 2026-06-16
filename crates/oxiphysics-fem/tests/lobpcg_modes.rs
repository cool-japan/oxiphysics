// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Integration tests for the LOBPCG shift-invert eigensolver.
//!
//! Ground truth is the 1-D Laplacian (`n = 100`, diagonal `2`, off-diagonal `-1`)
//! whose eigenvalues are known analytically: `λ_k = 4 sin²(k π / (2 (n+1)))`,
//! ascending for `k = 1, 2, …`.
//!
//! All tests are deterministic: the solver seeds its initial block with a fixed
//! seed and no global RNG is used.

use oxiphysics_fem::eigenvalue::{LobpcgConfig, lobpcg_solve};
use oxiphysics_fem::parallel_solver::CsrMatrix;

/// Build the 1-D Laplacian (`n × n`, tridiagonal, diag `2`, off-diag `-1`).
fn laplacian_1d(n: usize) -> CsrMatrix {
    let mut row_offsets = vec![0usize; n + 1];
    let mut col_indices: Vec<usize> = Vec::new();
    let mut values: Vec<f64> = Vec::new();
    for i in 0..n {
        if i > 0 {
            col_indices.push(i - 1);
            values.push(-1.0);
        }
        col_indices.push(i);
        values.push(2.0);
        if i + 1 < n {
            col_indices.push(i + 1);
            values.push(-1.0);
        }
        row_offsets[i + 1] = col_indices.len();
    }
    CsrMatrix {
        nrows: n,
        ncols: n,
        row_offsets,
        col_indices,
        values,
    }
}

/// Analytic eigenvalues of the 1-D Laplacian, ascending. `out[i]` is `λ_{i+1}`.
fn analytic_eigenvalues(n: usize) -> Vec<f64> {
    let denom = 2.0 * (n as f64 + 1.0);
    (1..=n)
        .map(|k| {
            let s = ((k as f64) * std::f64::consts::PI / denom).sin();
            4.0 * s * s
        })
        .collect()
}

/// HEADLINE GATE: first 20 modes within 0.5% of analytic (standard mode + AMG).
#[test]
fn plate_modes_within_half_percent() {
    let n = 100;
    let a = laplacian_1d(n);
    let analytic = analytic_eigenvalues(n);

    let cfg = LobpcgConfig {
        num_eigenvalues: 20,
        max_iter: 300,
        tol: 1e-6,
        shift: None,
        use_amg: true,
    };
    let res = lobpcg_solve(&a, None, &cfg).expect("LOBPCG standard solve returned an error");

    assert!(
        res.converged,
        "LOBPCG did not converge ({} iterations)",
        res.iterations
    );
    assert_eq!(res.eigenvalues.len(), 20, "expected 20 eigenpairs");

    for (i, &analytic_i) in analytic.iter().enumerate().take(20) {
        let rel = (res.eigenvalues[i] - analytic_i).abs() / analytic_i;
        assert!(
            rel < 0.005,
            "mode {i}: computed {} vs analytic {} (relative error {:.3e})",
            res.eigenvalues[i],
            analytic_i,
            rel
        );
    }

    assert_eq!(res.residual_norms.len(), 20, "expected 20 residual norms");
    for (i, &rn) in res.residual_norms.iter().enumerate() {
        assert!(
            rn.is_finite() && rn < 1e-3,
            "residual_norms[{i}] = {rn:.3e} not finite/below 1e-3"
        );
    }
}

/// The returned eigenvectors are orthonormal (B = I): ‖XᵀX − I‖_F < 1e-8.
#[test]
fn lobpcg_b_orthonormal() {
    let n = 100;
    let a = laplacian_1d(n);
    let k = 10;

    let cfg = LobpcgConfig {
        num_eigenvalues: k,
        max_iter: 300,
        tol: 1e-6,
        shift: None,
        use_amg: true,
    };
    let res = lobpcg_solve(&a, None, &cfg).expect("LOBPCG solve returned an error");
    assert_eq!(res.eigenvectors.len(), k, "expected {k} eigenvectors");

    let mut frob_sq = 0.0f64;
    for i in 0..k {
        for j in 0..k {
            let g: f64 = res.eigenvectors[i]
                .iter()
                .zip(res.eigenvectors[j].iter())
                .map(|(a_, b_)| a_ * b_)
                .sum();
            let target = if i == j { 1.0 } else { 0.0 };
            frob_sq += (g - target) * (g - target);
        }
    }
    let frob = frob_sq.sqrt();
    assert!(frob < 1e-8, "||XᵀX - I||_F = {frob:.3e}");
}

/// Shift-invert genuinely targets the interior: the 3 eigenpairs nearest σ = λ_10
/// land inside [λ_6, λ_14] and are NOT the three lowest modes.
#[test]
fn shift_invert_interior_targeting() {
    let n = 100;
    let a = laplacian_1d(n);
    let analytic = analytic_eigenvalues(n);
    let sigma = analytic[9]; // λ_10 (1-indexed 10th smallest)

    let cfg = LobpcgConfig {
        num_eigenvalues: 3,
        max_iter: 300,
        tol: 1e-6,
        shift: Some(sigma),
        use_amg: false,
    };
    let res = lobpcg_solve(&a, None, &cfg).expect("LOBPCG interior solve returned an error");
    assert_eq!(res.eigenvalues.len(), 3, "expected 3 interior eigenpairs");

    let lo = analytic[5]; // λ_6
    let hi = analytic[13]; // λ_14
    let lambda5 = analytic[4]; // λ_5

    for &ev in &res.eigenvalues {
        assert!(
            ev >= lo && ev <= hi,
            "eigenvalue {ev} is outside the interior window [{lo}, {hi}]; returned = {:?}",
            res.eigenvalues
        );
    }
    // Result is ascending, so the minimum is the first entry.
    let min_ev = res.eigenvalues[0];
    assert!(
        min_ev > lambda5,
        "minimum returned eigenvalue {min_ev} <= λ_5 {lambda5}: interior targeting failed \
         (looks like the lowest 3); returned = {:?}",
        res.eigenvalues
    );
}

/// AMG accelerates LOBPCG over Jacobi on the 1-D Laplacian: both converge and AMG
/// uses strictly fewer outer iterations.
#[test]
fn amg_preconditioner_efficacy() {
    let n = 100;
    let a = laplacian_1d(n);

    let cfg = |use_amg: bool| LobpcgConfig {
        num_eigenvalues: 5,
        max_iter: 300,
        tol: 1e-6,
        shift: None,
        use_amg,
    };

    let amg =
        lobpcg_solve(&a, None, &cfg(true)).expect("AMG-preconditioned solve returned an error");
    let jacobi =
        lobpcg_solve(&a, None, &cfg(false)).expect("Jacobi-preconditioned solve returned an error");

    assert!(amg.converged, "AMG run did not converge");
    assert!(
        jacobi.converged,
        "Jacobi run did not converge ({} iterations)",
        jacobi.iterations
    );
    // On the 1-D Laplacian a true multilevel V-cycle beats the (constant-diagonal)
    // Jacobi sweep, which behaves like no preconditioner — so AMG wins strictly.
    assert!(
        amg.iterations < jacobi.iterations,
        "expected AMG ({}) to use strictly fewer iterations than Jacobi ({})",
        amg.iterations,
        jacobi.iterations
    );
}
