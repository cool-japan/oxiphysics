// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Validates OxiPhysics FEM PCG+AMG against published algebraic multigrid benchmarks.
//!
//! Reference: Stuben 2001 "A review of algebraic multigrid"
//!            J. Comput. Appl. Math. 128:281-309.
//! Reference: Ruge & Stuben 1987 "Algebraic multigrid" in
//!            Multigrid Methods, SIAM (Frontiers in Applied Mathematics, Vol. 3).
//!
//! Key property validated: MESH-INDEPENDENCE — iteration count stays bounded as the grid
//! is refined, unlike classical iterative methods (PCG+Jacobi shows O(h⁻¹) growth).

use oxiphysics_fem::{
    parallel_solver::{CsrMatrix, ParallelPcgSolver, PcgWithPrecond},
    solvers::amg::{classical::AmgClassical, cycle::CycleKind, preconditioner::AmgPreconditioner},
};

// ── Convergence tolerance ─────────────────────────────────────────────────────

/// Relative residual tolerance used in all comparison tests.
const TOL: f64 = 1e-8;

/// Maximum PCG+AMG iterations allowed — well above the expected range so a
/// genuine non-convergence failure is distinguishable from an assertion failure.
const MAX_AMG_ITERS: usize = 200;

/// Maximum PCG+Jacobi iterations — generous so the solver actually converges
/// and the speedup ratio is meaningful.
const MAX_JAC_ITERS: usize = 2000;

// ── Matrix builders ───────────────────────────────────────────────────────────

/// Build a 3D finite-difference Poisson matrix on an n×n×n regular grid.
///
/// DOF index for node (i, j, k) = i*n*n + j*n + k.
/// 7-point stencil: diagonal = 6.0, six axis-aligned neighbours = -1.0.
/// Dirichlet boundary conditions are enforced by omitting off-grid couplings
/// (the boundary rows have diagonal > 6 as fewer neighbours are included,
/// which is the standard approach for finite-difference Poisson).
fn make_3d_poisson(n: usize) -> CsrMatrix {
    let ndofs = n * n * n;
    let mut row_offsets = vec![0usize; ndofs + 1];
    let mut col_indices: Vec<usize> = Vec::new();
    let mut values: Vec<f64> = Vec::new();

    for iz in 0..n {
        for iy in 0..n {
            for ix in 0..n {
                let dof = iz * n * n + iy * n + ix;

                // Collect neighbour entries sorted by global index.
                let mut entries: Vec<(usize, f64)> = Vec::with_capacity(7);

                // z-direction neighbours
                if iz > 0 {
                    entries.push((dof - n * n, -1.0));
                }
                // y-direction neighbours
                if iy > 0 {
                    entries.push((dof - n, -1.0));
                }
                // x-direction neighbours
                if ix > 0 {
                    entries.push((dof - 1, -1.0));
                }
                // diagonal = 6 always (standard 3D Poisson stencil with Dirichlet BCs;
                // missing off-grid couplings contribute 0 from the Dirichlet condition,
                // but the diagonal is defined by the underlying differential operator).
                entries.push((dof, 6.0));
                if ix + 1 < n {
                    entries.push((dof + 1, -1.0));
                }
                if iy + 1 < n {
                    entries.push((dof + n, -1.0));
                }
                if iz + 1 < n {
                    entries.push((dof + n * n, -1.0));
                }

                entries.sort_unstable_by_key(|&(c, _)| c);
                for (c, v) in entries {
                    col_indices.push(c);
                    values.push(v);
                }
                row_offsets[dof + 1] = col_indices.len();
            }
        }
    }

    CsrMatrix {
        nrows: ndofs,
        ncols: ndofs,
        row_offsets,
        col_indices,
        values,
    }
}

// ── Solver helpers ────────────────────────────────────────────────────────────

/// Run PCG+AMG on `(k, f)` and return `(iterations, converged)`.
///
/// Uses a fresh AMG hierarchy built from `k`.  The AMG preconditioner applies
/// one V-cycle per PCG iteration.
fn solve_pcg_amg(k: &CsrMatrix, f: &[f64]) -> (usize, bool) {
    let n = k.nrows;
    let amg_builder = AmgClassical::new();
    let hierarchy = amg_builder.build(k);

    let coarse_pcg = ParallelPcgSolver::new(500, 1e-10);
    let precond = AmgPreconditioner {
        hierarchy,
        cycle_kind: CycleKind::V,
        pcg: coarse_pcg,
    };

    let pcg = PcgWithPrecond::new(precond, MAX_AMG_ITERS, TOL);
    let mut x = vec![0.0f64; n];
    let stats = pcg.solve(k, f, &mut x);
    (stats.iterations, stats.converged)
}

/// Run PCG+Jacobi on `(k, f)` and return `(iterations, converged)`.
///
/// Uses the `ParallelPcgSolver` with diagonal (Jacobi) preconditioning.
fn solve_pcg_jacobi(k: &CsrMatrix, f: &[f64]) -> (usize, bool) {
    let n = k.nrows;
    let pcg = ParallelPcgSolver::new(MAX_JAC_ITERS, TOL);
    let mut x = vec![0.0f64; n];
    let stats = pcg.solve(k, f, &mut x);
    (stats.iterations, stats.converged)
}

// ── Tests ─────────────────────────────────────────────────────────────────────

/// PCG+AMG on an 8³ 3D Poisson system must converge within the published
/// iteration range.
///
/// Stuben 2001 (Table 3) reports ~7–10 iterations for standard RS-AMG on 3D
/// Poisson at tolerance 1e-8.  We allow a 2× factor to account for
/// implementation-specific differences, giving an upper bound of 24 iterations.
#[test]
fn test_amg_converges_within_published_iteration_range() {
    let n = 8;
    let k = make_3d_poisson(n);
    let f = vec![1.0f64; k.nrows];

    let (iters, converged) = solve_pcg_amg(&k, &f);

    assert!(
        converged,
        "PCG+AMG on {}³ Poisson did not converge within {MAX_AMG_ITERS} iters (got {iters} iters)",
        n
    );
    assert!(
        iters >= 1,
        "PCG+AMG on {}³ Poisson: solver did not advance (iters=0)",
        n
    );
    assert!(
        iters <= 24,
        "PCG+AMG on {}³ Poisson: expected ≤24 iters (Stuben 2001 reference: 7–10), got {iters}",
        n
    );

    eprintln!("AMG {}³ Poisson: {iters} iterations (reference: 7–10)", n);
}

/// PCG+AMG must exhibit mesh-independent convergence: the iteration count must
/// not grow significantly as the grid is refined from 8³ to 32³.
///
/// This is the defining property of algebraic multigrid (Stuben 2001, Sec. 3).
/// Classical PCG+Jacobi shows O(h⁻¹) growth; true AMG should show O(1) or at
/// most O(log(h⁻¹)) behaviour.  We allow at most 6 additional iterations across
/// a 4× refinement in each dimension (512→32 768 DOF, factor 64 in problem size).
#[test]
fn test_amg_mesh_independence() {
    let f_8 = vec![1.0f64; 8 * 8 * 8];
    let f_16 = vec![1.0f64; 16 * 16 * 16];
    let f_32 = vec![1.0f64; 32 * 32 * 32];

    let (iters_8, conv_8) = solve_pcg_amg(&make_3d_poisson(8), &f_8);
    let (iters_16, conv_16) = solve_pcg_amg(&make_3d_poisson(16), &f_16);
    let (iters_32, conv_32) = solve_pcg_amg(&make_3d_poisson(32), &f_32);

    assert!(
        conv_8,
        "PCG+AMG on 8³ Poisson did not converge (iters={iters_8})"
    );
    assert!(
        conv_16,
        "PCG+AMG on 16³ Poisson did not converge (iters={iters_16})"
    );
    assert!(
        conv_32,
        "PCG+AMG on 32³ Poisson did not converge (iters={iters_32})"
    );

    // Mesh-independence criterion: growth from coarsest to finest ≤ 6 iters.
    // (Stuben 2001 reports growth ≤ 4 for ideal RS-AMG; we allow ≤ 6.)
    let growth = iters_32.saturating_sub(iters_8);
    assert!(
        growth <= 6,
        "AMG mesh-independence FAILED: iteration count grew by {growth} from 8³({iters_8}) \
         to 32³({iters_32}). Reference (Stuben 2001): growth should be ≤4 for true AMG. \
         16³ count was {iters_16}."
    );

    eprintln!(
        "AMG mesh-independence: 8³={iters_8}, 16³={iters_16}, 32³={iters_32} (growth={growth} iters)"
    );
}

/// PCG+AMG must use significantly fewer iterations than PCG+Jacobi.
///
/// The existing `perf_bench::amg_poisson` benchmark confirms PCG+Jacobi needs
/// ≥200 iterations on 128² 2D Poisson at tol 1e-8.  For the 3D 32³ problem we
/// require AMG to use at most 1/5 of the Jacobi iteration count (≥5× speedup).
#[test]
fn test_amg_vs_jacobi_speedup() {
    let n = 32;
    let k = make_3d_poisson(n);
    let f = vec![1.0f64; k.nrows];

    let (iters_amg, conv_amg) = solve_pcg_amg(&k, &f);
    let (iters_jac, conv_jac) = solve_pcg_jacobi(&k, &f);

    assert!(
        conv_amg,
        "PCG+AMG on {}³ Poisson did not converge (iters={iters_amg})",
        n
    );
    assert!(
        conv_jac,
        "PCG+Jacobi on {}³ Poisson did not converge within {MAX_JAC_ITERS} iters \
         (iters={iters_jac}); increase MAX_JAC_ITERS or reduce grid size",
        n
    );

    // Require AMG to use at most 1/5 of the Jacobi iteration count.
    let speedup = iters_jac as f64 / iters_amg.max(1) as f64;
    assert!(
        speedup >= 5.0,
        "AMG speedup over Jacobi: expected ≥5× on {}³ Poisson, got {:.1}× \
         ({iters_amg} AMG iters vs {iters_jac} Jacobi iters)",
        n,
        speedup
    );

    eprintln!(
        "AMG vs Jacobi on {}³ Poisson: {iters_amg} vs {iters_jac} iters = {:.1}× speedup",
        n, speedup
    );
}
