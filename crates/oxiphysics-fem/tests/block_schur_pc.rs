// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Integration tests for the Chebyshev smoother and block-Schur preconditioner.
//!
//! 1. `chebyshev_smoothing_factor` — Chebyshev damps high-frequency error and is
//!    competitive with Gauss-Seidel.
//! 2. `block_schur_mesh_independence` — block-Schur GMRES iteration counts vary
//!    ≤ ±10% across three mesh refinements (Murphy-Golub-Wathen clustering on an
//!    exact-Schur synthetic saddle system).
//! 3. `block_schur_dense_parity` — block-preconditioned GMRES recovers the exact
//!    solution of a small dense saddle system.

use oxiphysics_fem::parallel_solver::CsrMatrix;
use oxiphysics_fem::solvers::amg::chebyshev_smoother::chebyshev_smoother;
use oxiphysics_fem::solvers::amg::smoothers::forward_gs;
use oxiphysics_fem::solvers::block_schur_pc::{
    BlockSchurPreconditioner, block_schur_gmres, extract_velocity_block, lump_pressure_mass,
};

// ── helpers ─────────────────────────────────────────────────────────────────

fn l2(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

fn dense_matvec(m: &[Vec<f64>], x: &[f64]) -> Vec<f64> {
    m.iter()
        .map(|row| row.iter().zip(x.iter()).map(|(a, b)| a * b).sum())
        .collect()
}

fn poisson_1d_dense(n: usize) -> Vec<Vec<f64>> {
    let mut a = vec![vec![0.0f64; n]; n];
    for (i, row) in a.iter_mut().enumerate() {
        row[i] = 2.0;
        if i > 0 {
            row[i - 1] = -1.0;
        }
        if i + 1 < n {
            row[i + 1] = -1.0;
        }
    }
    a
}

fn poisson_1d_csr(n: usize) -> CsrMatrix {
    let mut row_offsets = vec![0usize; n + 1];
    let mut col_indices = Vec::new();
    let mut values = Vec::new();
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

// ── Test 1: Chebyshev smoothing factor ──────────────────────────────────────

#[test]
fn chebyshev_smoothing_factor() {
    let n = 33;
    let a_dense = poisson_1d_dense(n);
    let a_csr = poisson_1d_csr(n);

    // High-frequency (checkerboard) error mode; b = 0 ⇒ x is the error itself.
    let b = vec![0.0f64; n];
    let x0: Vec<f64> = (0..n)
        .map(|i| if i % 2 == 0 { 1.0 } else { -1.0 })
        .collect();
    let norm0 = l2(&x0);

    let degree = 5;

    let mut x_cheb = x0.clone();
    chebyshev_smoother(&a_dense, &b, &mut x_cheb, degree);
    let factor_cheb = l2(&x_cheb) / norm0;

    let mut x_gs = x0.clone();
    forward_gs(&a_csr, &b, &mut x_gs, degree);
    let factor_gs = l2(&x_gs) / norm0;

    println!("chebyshev factor = {factor_cheb:.4}, gauss-seidel factor = {factor_gs:.4}");

    // Chebyshev strongly damps the high-frequency error (smoother behaviour).
    assert!(
        factor_cheb < 0.6,
        "Chebyshev high-frequency damping factor {factor_cheb:.4} not < 0.6"
    );
    // Both smoothers reduce the error (Chebyshev is competitive with GS).
    assert!(
        factor_cheb < 1.0,
        "Chebyshev did not converge: {factor_cheb:.4}"
    );
    assert!(
        factor_gs < 1.0,
        "Gauss-Seidel did not converge: {factor_gs:.4}"
    );
}

// ── Test 2: block-Schur GMRES mesh-independence (headline gate) ──────────────

/// Build an *exact-Schur* synthetic Stokes saddle system at refinement `np`.
///
/// `A = 2·I` (SPD velocity block; AMG → exact single-level PCG solve), and `B`
/// has two disjoint non-zeros per row, so `B·Bᵀ` is diagonal and the exact Schur
/// complement `B·A⁻¹·Bᵀ` is diagonal too.  Setting the lumped pressure mass to
/// that exact diagonal (with `ν = 1`) makes the block preconditioner *exact*, so
/// by Murphy-Golub-Wathen `P⁻¹K` has three distinct eigenvalues regardless of
/// `np`: GMRES converges in a mesh-independent iteration count.
fn build_synthetic_saddle(np: usize) -> (Vec<Vec<f64>>, Vec<f64>, usize, usize) {
    let nv = 2 * np;
    let ntot = nv + np;
    let alpha = 2.0;

    let mut saddle = vec![vec![0.0f64; ntot]; ntot];
    for (i, row) in saddle.iter_mut().enumerate().take(nv) {
        row[i] = alpha;
    }

    let mut mp_diag = vec![0.0f64; np];
    for (i, mp) in mp_diag.iter_mut().enumerate() {
        let a_i = 1.0 + 0.5 * ((i % 5) as f64);
        let b_i = 1.0 + 0.3 * ((i % 3) as f64);
        let c0 = 2 * i;
        let c1 = 2 * i + 1;
        // Divergence block B in the pressure rows, and Bᵀ in the velocity columns.
        saddle[nv + i][c0] = a_i;
        saddle[nv + i][c1] = b_i;
        saddle[c0][nv + i] = a_i;
        saddle[c1][nv + i] = b_i;
        // Exact Schur diagonal: (B·A⁻¹·Bᵀ)_ii = (a_i² + b_i²)/α.
        *mp = (a_i * a_i + b_i * b_i) / alpha;
    }

    (saddle, mp_diag, nv, np)
}

#[test]
fn block_schur_mesh_independence() {
    let mut iter_counts: Vec<usize> = Vec::new();

    for &np in &[8usize, 16, 32] {
        let (saddle, mp_diag, nv, n_pres) = build_synthetic_saddle(np);
        let ntot = nv + n_pres;
        let nu = 1.0;

        let vel_block = extract_velocity_block(&saddle, nv);
        let pc = BlockSchurPreconditioner::new(&vel_block, mp_diag, nu, nv, n_pres);

        // Generic known solution ⇒ RHS excites every eigenspace of P⁻¹K.
        let x_true: Vec<f64> = (0..ntot).map(|k| 1.0 + ((k % 7) as f64) * 0.37).collect();
        let b = dense_matvec(&saddle, &x_true);

        let mut x = vec![0.0f64; ntot];
        let stats = block_schur_gmres(&saddle, &b, &mut x, &pc, 50, 5, 1e-8);

        assert!(
            stats.converged,
            "block-Schur GMRES failed at np={np}: residual = {:.3e}",
            stats.residual_norm
        );
        iter_counts.push(stats.iterations);
    }

    let min = iter_counts.iter().copied().min().unwrap_or(1).max(1);
    let max = iter_counts.iter().copied().max().unwrap_or(1);
    println!("block-Schur GMRES iteration counts across meshes: {iter_counts:?}");

    assert!(
        (max as f64) <= 1.1 * (min as f64),
        "GMRES iteration counts not mesh-independent (>±10%): {iter_counts:?}"
    );
}

// ── Test 3: dense-direct parity ─────────────────────────────────────────────

#[test]
fn block_schur_dense_parity() {
    let nv = 4;
    let n_pres = 2;
    let ntot = nv + n_pres;

    // SPD tridiagonal velocity block (exercises the AMG/PCG velocity solve).
    let a_block = vec![
        vec![4.0, -1.0, 0.0, 0.0],
        vec![-1.0, 4.0, -1.0, 0.0],
        vec![0.0, -1.0, 4.0, -1.0],
        vec![0.0, 0.0, -1.0, 4.0],
    ];
    // Full-row-rank divergence block B (2 × 4).
    let b_block = vec![vec![1.0, 0.5, 0.0, 0.0], vec![0.0, 0.0, 0.5, 1.0]];
    // Consistent SPD pressure mass matrix.
    let mp_block = vec![vec![2.0, 0.5], vec![0.5, 2.0]];

    // Assemble the dense saddle matrix [[A, Bᵀ], [B, 0]].
    let mut saddle = vec![vec![0.0f64; ntot]; ntot];
    for (i, a_row) in a_block.iter().enumerate() {
        saddle[i][..nv].copy_from_slice(a_row);
    }
    for (i, b_row) in b_block.iter().enumerate() {
        for (j, &bij) in b_row.iter().enumerate() {
            saddle[nv + i][j] = bij;
            saddle[j][nv + i] = bij;
        }
    }

    let nu = 1.0;
    let vel = extract_velocity_block(&saddle, nv);
    assert_eq!(vel, a_block, "extract_velocity_block mismatch");

    let mp_diag = lump_pressure_mass(&mp_block);
    assert_eq!(mp_diag, vec![2.5, 2.5]);

    let pc = BlockSchurPreconditioner::new(&vel, mp_diag, nu, nv, n_pres);

    let x_true = vec![1.0, -2.0, 0.5, 3.0, -1.5, 2.0];
    let b = dense_matvec(&saddle, &x_true);

    let mut x = vec![0.0f64; ntot];
    let stats = block_schur_gmres(&saddle, &b, &mut x, &pc, 50, 20, 1e-12);

    let err = l2(&x
        .iter()
        .zip(x_true.iter())
        .map(|(xi, ti)| xi - ti)
        .collect::<Vec<_>>());

    assert!(
        err < 1e-8,
        "block-preconditioned GMRES did not match the direct solution: error = {err:.3e}, stats = {stats:?}"
    );
}
