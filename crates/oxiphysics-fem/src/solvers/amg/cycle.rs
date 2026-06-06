// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! V-cycle and W-cycle for AMG hierarchies.

use crate::parallel_solver::{CsrMatrix, ParallelPcgSolver, PcgStats};
use crate::solvers::amg::smoothers::symmetric_gs;

/// One level in the AMG hierarchy.
#[derive(Debug, Clone)]
pub struct AmgLevel {
    /// System matrix at this level.
    pub a: CsrMatrix,
    /// Prolongation operator (coarse → fine).
    pub p: CsrMatrix,
    /// Restriction operator P^T (fine → coarse).
    pub pt: CsrMatrix,
    /// Diagonal entries of `a` (for reference; GS smoother uses a directly).
    pub diag: Vec<f64>,
}

/// The complete AMG level hierarchy.
#[derive(Debug, Clone)]
pub struct AmgHierarchy {
    /// Levels from finest (index 0) to coarsest (last).
    pub levels: Vec<AmgLevel>,
    /// Problem size below which we call the coarse solver directly.
    pub coarse_cutoff: usize,
}

/// V-cycle or W-cycle selector.
#[derive(Debug, Clone, Copy)]
pub enum CycleKind {
    /// Standard V-cycle (one coarse-grid correction per level).
    V,
    /// W-cycle (two coarse-grid corrections per level, more aggressive).
    W,
}

/// Run one V-cycle at the given level.
///
/// - Pre-smooth with symmetric GS.
/// - Restrict residual to coarse level.
/// - Recursively correct.
/// - Prolongate correction and add.
/// - Post-smooth with symmetric GS.
pub fn v_cycle(
    hier: &AmgHierarchy,
    level: usize,
    b: &[f64],
    x: &mut [f64],
    pcg: &ParallelPcgSolver,
) {
    let last_level = hier.levels.len() - 1;
    let lev = &hier.levels[level];

    // Coarse solve: use PCG directly
    if level == last_level || lev.a.nrows <= hier.coarse_cutoff {
        pcg.solve(&lev.a, b, x);
        return;
    }

    // Pre-smooth
    symmetric_gs(&lev.a, b, x, 2);

    // Compute residual r = b - A*x
    let n = lev.a.nrows;
    let mut ax = vec![0.0f64; n];
    lev.a.spmv(x, &mut ax);
    let r: Vec<f64> = b.iter().zip(ax.iter()).map(|(bi, ai)| bi - ai).collect();

    // Restrict residual: r_c = P^T * r
    let n_coarse = hier.levels[level + 1].a.nrows;
    let mut r_c = vec![0.0f64; n_coarse];
    lev.pt.spmv(&r, &mut r_c);

    // Coarse correction
    let mut x_c = vec![0.0f64; n_coarse];
    v_cycle(hier, level + 1, &r_c, &mut x_c, pcg);

    // Prolongate and add: x += P * x_c
    let mut p_xc = vec![0.0f64; n];
    lev.p.spmv(&x_c, &mut p_xc);
    for i in 0..n {
        x[i] += p_xc[i];
    }

    // Post-smooth
    symmetric_gs(&lev.a, b, x, 2);
}

/// Run one W-cycle at the given level.
///
/// Same as V-cycle but applies the coarse correction twice:
/// - First correction on the restricted residual.
/// - Second correction on the residual of the first correction.
pub fn w_cycle(
    hier: &AmgHierarchy,
    level: usize,
    b: &[f64],
    x: &mut [f64],
    pcg: &ParallelPcgSolver,
) {
    let last_level = hier.levels.len() - 1;
    let lev = &hier.levels[level];

    // Coarse solve
    if level == last_level || lev.a.nrows <= hier.coarse_cutoff {
        pcg.solve(&lev.a, b, x);
        return;
    }

    // Pre-smooth
    symmetric_gs(&lev.a, b, x, 2);

    // Compute residual
    let n = lev.a.nrows;
    let mut ax = vec![0.0f64; n];
    lev.a.spmv(x, &mut ax);
    let r: Vec<f64> = b.iter().zip(ax.iter()).map(|(bi, ai)| bi - ai).collect();

    // Restrict residual
    let n_coarse = hier.levels[level + 1].a.nrows;
    let mut r_c = vec![0.0f64; n_coarse];
    lev.pt.spmv(&r, &mut r_c);

    // First coarse correction
    let mut x_c = vec![0.0f64; n_coarse];
    w_cycle(hier, level + 1, &r_c, &mut x_c, pcg);

    // Second coarse correction: compute residual of x_c, then correct
    let mut a_c_xc = vec![0.0f64; n_coarse];
    hier.levels[level + 1].a.spmv(&x_c, &mut a_c_xc);
    let r_c2: Vec<f64> = r_c
        .iter()
        .zip(a_c_xc.iter())
        .map(|(ri, ai)| ri - ai)
        .collect();
    let mut x_c2 = vec![0.0f64; n_coarse];
    w_cycle(hier, level + 1, &r_c2, &mut x_c2, pcg);
    for i in 0..n_coarse {
        x_c[i] += x_c2[i];
    }

    // Prolongate and add
    let mut p_xc = vec![0.0f64; n];
    lev.p.spmv(&x_c, &mut p_xc);
    for i in 0..n {
        x[i] += p_xc[i];
    }

    // Post-smooth
    symmetric_gs(&lev.a, b, x, 2);
}

/// Solve using the AMG hierarchy with V-cycles or W-cycles.
///
/// Runs up to `max_cycles` of the chosen cycle type, stopping when the
/// relative residual `||r|| / ||b||` drops below `tol`.
pub fn amg_solve(
    hier: &AmgHierarchy,
    b: &[f64],
    x: &mut [f64],
    kind: CycleKind,
    max_cycles: usize,
    tol: f64,
) -> PcgStats {
    let n = hier.levels[0].a.nrows;
    let pcg = ParallelPcgSolver::new(500, 1e-10);

    let b_norm = b.iter().map(|bi| bi * bi).sum::<f64>().sqrt().max(1e-300);
    let mut res_norm;

    for cycle in 0..max_cycles {
        // Compute current residual
        let mut ax = vec![0.0f64; n];
        hier.levels[0].a.spmv(x, &mut ax);
        let r: Vec<f64> = b.iter().zip(ax.iter()).map(|(bi, ai)| bi - ai).collect();
        res_norm = r.iter().map(|ri| ri * ri).sum::<f64>().sqrt();

        if res_norm / b_norm < tol {
            return PcgStats {
                iterations: cycle,
                residual_norm: res_norm,
                converged: true,
            };
        }

        match kind {
            CycleKind::V => v_cycle(hier, 0, b, x, &pcg),
            CycleKind::W => w_cycle(hier, 0, b, x, &pcg),
        }
    }

    // Final residual after all cycles
    let mut ax = vec![0.0f64; n];
    hier.levels[0].a.spmv(x, &mut ax);
    let r: Vec<f64> = b.iter().zip(ax.iter()).map(|(bi, ai)| bi - ai).collect();
    res_norm = r.iter().map(|ri| ri * ri).sum::<f64>().sqrt();

    PcgStats {
        iterations: max_cycles,
        residual_norm: res_norm,
        converged: res_norm / b_norm < tol,
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a 1D Poisson tridiagonal matrix.
    fn make_1d_poisson(n: usize) -> CsrMatrix {
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

    /// Build a trivial 2-level hierarchy manually for testing v_cycle.
    fn make_trivial_2level_hierarchy(n_fine: usize) -> AmgHierarchy {
        // Fine level: 1D Poisson of size n_fine
        let a_fine = make_1d_poisson(n_fine);

        // C-points: even-indexed nodes; F-points: odd-indexed
        // n_coarse = ceil(n_fine / 2)
        let n_coarse = n_fine.div_ceil(2);

        // Build P: C-point rows get weight 1.0; F-point rows interpolate from neighbors
        let mut row_offsets = vec![0usize; n_fine + 1];
        let mut col_indices = Vec::new();
        let mut values = Vec::new();

        for i in 0..n_fine {
            let c_idx = i / 2; // coarse index for even-indexed nodes
            if i % 2 == 0 {
                // C-point: direct injection
                col_indices.push(c_idx);
                values.push(1.0);
                row_offsets[i + 1] = col_indices.len();
            } else {
                // F-point: average of neighbors
                let left_c = i / 2; // left C-point index (even, left)
                let right_c = i.div_ceil(2); // right C-point index
                if right_c < n_coarse {
                    col_indices.push(left_c);
                    values.push(0.5);
                    col_indices.push(right_c);
                    values.push(0.5);
                } else {
                    col_indices.push(left_c);
                    values.push(1.0);
                }
                row_offsets[i + 1] = col_indices.len();
            }
        }

        let p = CsrMatrix {
            nrows: n_fine,
            ncols: n_coarse,
            row_offsets,
            col_indices,
            values,
        };

        use crate::solvers::amg::galerkin::{csr_transpose, galerkin_coarse};
        let pt = csr_transpose(&p);
        let a_coarse = galerkin_coarse(&a_fine, &p);

        let diag_fine: Vec<f64> = (0..n_fine)
            .map(|i| {
                let rs = a_fine.row_offsets[i];
                let re = a_fine.row_offsets[i + 1];
                let mut d = 1.0;
                for k in rs..re {
                    if a_fine.col_indices[k] == i {
                        d = a_fine.values[k];
                        break;
                    }
                }
                d
            })
            .collect();

        let diag_coarse: Vec<f64> = (0..n_coarse)
            .map(|i| {
                let rs = a_coarse.row_offsets[i];
                let re = a_coarse.row_offsets[i + 1];
                let mut d = 1.0;
                for k in rs..re {
                    if a_coarse.col_indices[k] == i {
                        d = a_coarse.values[k];
                        break;
                    }
                }
                d
            })
            .collect();

        AmgHierarchy {
            levels: vec![
                AmgLevel {
                    a: a_fine,
                    p,
                    pt,
                    diag: diag_fine,
                },
                AmgLevel {
                    a: a_coarse,
                    p: CsrMatrix::identity(0),
                    pt: CsrMatrix::identity(0),
                    diag: diag_coarse,
                },
            ],
            coarse_cutoff: 4,
        }
    }

    #[test]
    fn test_vcycle_1d_poisson_reduces() {
        // Test that one V-cycle on 1D Poisson reduces the residual
        let n = 32;
        let hier = make_trivial_2level_hierarchy(n);
        let a = &hier.levels[0].a;

        let b: Vec<f64> = (0..n).map(|i| (i + 1) as f64).collect();
        let mut x = vec![0.0f64; n];

        // Compute initial residual
        let r0_norm: f64 = b.iter().map(|bi| bi * bi).sum::<f64>().sqrt();

        let pcg = ParallelPcgSolver::new(500, 1e-10);
        v_cycle(&hier, 0, &b, &mut x, &pcg);

        let mut ax = vec![0.0f64; n];
        a.spmv(&x, &mut ax);
        let r1_norm: f64 = b
            .iter()
            .zip(ax.iter())
            .map(|(bi, ai)| (bi - ai).powi(2))
            .sum::<f64>()
            .sqrt();

        let ratio = r1_norm / r0_norm;
        assert!(
            ratio <= 0.5,
            "V-cycle residual reduction ratio {ratio:.4} exceeds 0.5"
        );
    }
}
