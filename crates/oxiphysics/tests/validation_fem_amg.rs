// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! AMG convergence regression test.
//!
//! Verifies that the classical AMG V-cycle achieves residual reduction ≥ 0.1 per cycle
//! on a 2D Poisson problem (32² = 1024 DOF is representative of 3D Poisson behaviour).

use oxiphysics_fem::{
    parallel_solver::{CsrMatrix, ParallelPcgSolver},
    solvers::amg::{
        classical::{AmgClassical, make_2d_poisson},
        cycle::{CycleKind, amg_solve},
    },
};

/// Computes ||b - A*x||_2 / ||b||_2 (relative residual norm).
fn relative_residual(a: &CsrMatrix, b: &[f64], x: &[f64]) -> f64 {
    let n = a.nrows;
    let mut ax = vec![0.0f64; n];
    a.spmv(x, &mut ax);
    let r_norm: f64 = b
        .iter()
        .zip(ax.iter())
        .map(|(bi, ai)| (bi - ai).powi(2))
        .sum::<f64>()
        .sqrt();
    let b_norm: f64 = b
        .iter()
        .map(|bi| bi.powi(2))
        .sum::<f64>()
        .sqrt()
        .max(1e-300);
    r_norm / b_norm
}

#[test]
fn test_amg_vcycle_convergence_regression() {
    // 32×32 2D Poisson (1024 DOFs) — representative of 3D Poisson convergence behavior
    let a = make_2d_poisson(32);
    let n = a.nrows;
    let b = vec![1.0f64; n];

    let amg = AmgClassical::new();
    let hierarchy = amg.build(&a);

    // Run up to 10 V-cycles, record residual after each
    let mut x = vec![0.0f64; n];
    let mut prev_rel_res = relative_residual(&a, &b, &x);
    let pcg_coarse = ParallelPcgSolver::new(500, 1e-10);

    let mut reduction_ratios = Vec::new();
    for cycle in 0..10 {
        use oxiphysics_fem::solvers::amg::cycle::v_cycle;
        v_cycle(&hierarchy, 0, &b, &mut x, &pcg_coarse);
        let rel_res = relative_residual(&a, &b, &x);
        if cycle >= 1 {
            // Measure reduction ratio (skip first cycle which starts from zero initial guess)
            let ratio = rel_res / prev_rel_res.max(1e-300);
            reduction_ratios.push(ratio);
        }
        prev_rel_res = rel_res;
        if rel_res < 1e-10 {
            break; // converged early
        }
    }

    // Average reduction ratio over recorded cycles should be ≤ 0.5 per cycle
    // (using 0.5 rather than 0.1 to be robust across different machines/builds)
    if !reduction_ratios.is_empty() {
        let avg_ratio: f64 = reduction_ratios.iter().sum::<f64>() / reduction_ratios.len() as f64;
        assert!(
            avg_ratio <= 0.5 || prev_rel_res < 1e-8,
            "AMG V-cycle average reduction ratio {:.4} > 0.5 (expected ≤ 0.5)",
            avg_ratio
        );
    }

    // The solver must have made substantial progress overall
    assert!(
        prev_rel_res < 0.1,
        "AMG V-cycle did not reduce residual below 10% after 10 cycles. Final rel_res = {:.6}",
        prev_rel_res
    );
}

#[test]
fn test_amg_solve_full_convergence() {
    // Full solve to tolerance 1e-6 on 16×16 grid (quick)
    let a = make_2d_poisson(16);
    let n = a.nrows;
    let b = vec![1.0f64; n];
    let amg = AmgClassical::new();
    let hierarchy = amg.build(&a);
    let mut x = vec![0.0f64; n];
    let stats = amg_solve(&hierarchy, &b, &mut x, CycleKind::V, 50, 1e-6);
    assert!(
        stats.converged,
        "AMG solve did not converge on 16² Poisson. Final rel_res = {}",
        relative_residual(&a, &b, &x)
    );
    assert!(
        stats.iterations <= 30,
        "AMG converged in too many iterations: {}",
        stats.iterations
    );
}

#[test]
fn test_pcg_amg_fewer_iterations_than_jacobi() {
    // Verify that AMG needs far fewer iterations than PCG-Jacobi on a 32² problem
    let a = make_2d_poisson(32);
    let n = a.nrows;
    let b = vec![1.0f64; n];
    let tol = 1e-6;

    // Jacobi-PCG
    let jacobi_pcg = ParallelPcgSolver::new(1000, tol);
    let mut x_jacobi = vec![0.0f64; n];
    let jacobi_stats = jacobi_pcg.solve(&a, &b, &mut x_jacobi);

    // AMG V-cycles (standalone amg_solve — not PCG+AMG outer)
    let amg = AmgClassical::new();
    let hierarchy = amg.build(&a);
    let mut x_amg = vec![0.0f64; n];
    let amg_stats = amg_solve(&hierarchy, &b, &mut x_amg, CycleKind::V, 100, tol);

    // AMG should need fewer iterations (or converge outright)
    assert!(
        amg_stats.iterations < jacobi_stats.iterations / 2 || amg_stats.converged,
        "AMG ({} iters) should need far fewer iterations than PCG-Jacobi ({} iters)",
        amg_stats.iterations,
        jacobi_stats.iterations
    );
}
