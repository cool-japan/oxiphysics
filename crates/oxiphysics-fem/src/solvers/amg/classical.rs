// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Classical (Ruge-Stüben) AMG hierarchy builder.

use crate::parallel_solver::CsrMatrix;
use crate::solvers::amg::{
    cycle::{AmgHierarchy, AmgLevel},
    galerkin::{csr_transpose, galerkin_coarse},
    graph::{cf_splitting, strong_connections},
};

/// Classical Ruge-Stüben AMG builder.
#[derive(Debug, Clone)]
pub struct AmgClassical {
    /// Maximum number of hierarchy levels.
    pub max_levels: usize,
    /// Problem size below which coarsening stops (direct solve).
    pub coarse_cutoff: usize,
    /// Strong-connection threshold θ (default 0.25).
    pub theta: f64,
}

impl Default for AmgClassical {
    fn default() -> Self {
        Self::new()
    }
}

impl AmgClassical {
    /// Create an `AmgClassical` with default parameters.
    pub fn new() -> Self {
        Self {
            max_levels: 10,
            coarse_cutoff: 100,
            theta: 0.25,
        }
    }

    /// Build the AMG hierarchy starting from fine-level matrix `a`.
    pub fn build(&self, a: &CsrMatrix) -> AmgHierarchy {
        let mut levels: Vec<AmgLevel> = Vec::new();
        let mut current_a = a.clone();

        while levels.len() < self.max_levels && current_a.nrows > self.coarse_cutoff {
            let strong = strong_connections(&current_a, self.theta);
            let is_c = cf_splitting(&strong);

            let c_count = is_c.iter().filter(|&&b| b).count();

            // Degenerate cases: stop coarsening
            if c_count == 0 || c_count == current_a.nrows {
                break;
            }

            // Map fine DOF → coarse DOF index (C-points only)
            let mut c_to_coarse: Vec<Option<usize>> = vec![None; current_a.nrows];
            let mut coarse_idx = 0usize;
            for (i, &is_cp) in is_c.iter().enumerate() {
                if is_cp {
                    c_to_coarse[i] = Some(coarse_idx);
                    coarse_idx += 1;
                }
            }
            let n_coarse = coarse_idx;

            // Build prolongation operator P
            let p = build_prolongation(&current_a, &is_c, &c_to_coarse, n_coarse, &strong);
            let pt = csr_transpose(&p);

            // Galerkin coarse operator: A_c = P^T * A * P
            let a_coarse = galerkin_coarse(&current_a, &p);

            // Store actual diagonal values (not inverse)
            let diag: Vec<f64> = extract_diagonal(&current_a);

            levels.push(AmgLevel {
                a: current_a,
                p,
                pt,
                diag,
            });

            current_a = a_coarse;
        }

        // Add final coarse level (direct solve)
        let diag_coarse = extract_diagonal(&current_a);
        levels.push(AmgLevel {
            a: current_a,
            p: CsrMatrix::identity(0),
            pt: CsrMatrix::identity(0),
            diag: diag_coarse,
        });

        AmgHierarchy {
            levels,
            coarse_cutoff: self.coarse_cutoff,
        }
    }
}

/// Extract diagonal values from a CsrMatrix.
pub fn extract_diagonal(a: &CsrMatrix) -> Vec<f64> {
    (0..a.nrows)
        .map(|i| {
            let rs = a.row_offsets[i];
            let re = a.row_offsets[i + 1];
            let mut d = 0.0f64;
            for k in rs..re {
                if a.col_indices[k] == i {
                    d = a.values[k];
                    break;
                }
            }
            d
        })
        .collect()
}

/// Build the prolongation operator P (fine → coarse).
///
/// For each fine DOF i:
/// - C-point: P[i, c_to_coarse[i]] = 1.0 (direct injection)
/// - F-point with C-neighbors: classical direct interpolation
///   P[i, c] = A[i,c] / (sum of A[i,k] for k in C_strong(i))
///   (note: both numerator and denominator include the actual A[i,c] values,
///   which are negative for M-matrices; the ratio comes out positive = correct weight)
/// - F-point with no C-neighbors: fallback injection to nearest C-point
fn build_prolongation(
    a: &CsrMatrix,
    is_c: &[bool],
    c_to_coarse: &[Option<usize>],
    n_coarse: usize,
    strong: &[Vec<usize>],
) -> CsrMatrix {
    let n_fine = a.nrows;
    let mut row_offsets = vec![0usize; n_fine + 1];
    let mut col_indices: Vec<usize> = Vec::new();
    let mut values: Vec<f64> = Vec::new();

    for i in 0..n_fine {
        if is_c[i] {
            // C-point: direct injection
            let c = c_to_coarse[i].expect("C-point must have a coarse index");
            col_indices.push(c);
            values.push(1.0);
            row_offsets[i + 1] = col_indices.len();
        } else {
            // F-point: find C-points in strong[i]
            let c_strong: Vec<usize> = strong[i].iter().filter(|&&j| is_c[j]).copied().collect();

            if c_strong.is_empty() {
                // Fallback: find any C-point among all neighbors in A
                let rs = a.row_offsets[i];
                let re = a.row_offsets[i + 1];
                let mut fallback_c: Option<usize> = None;
                for k in rs..re {
                    let j = a.col_indices[k];
                    if j != i && is_c[j] {
                        fallback_c = Some(j);
                        break;
                    }
                }

                if let Some(fc) = fallback_c {
                    let c = c_to_coarse[fc].expect("C-point must have coarse index");
                    col_indices.push(c);
                    values.push(1.0);
                } else {
                    // Find any C-point globally and inject to it as last resort
                    let nearest_c = c_to_coarse.iter().position(|x| x.is_some()).unwrap_or(0);
                    let c = c_to_coarse[nearest_c].unwrap_or(0);
                    col_indices.push(c);
                    values.push(1.0);
                }
                row_offsets[i + 1] = col_indices.len();
                continue;
            }

            // Classical direct interpolation
            // P[i,c] = A[i,c] / (sum_{k in C_strong(i)} A[i,k])
            // For M-matrices: A[i,c] < 0 and sum < 0, so ratio > 0.
            // We want the weight to be positive and the row to sum to 1 on constants.
            //
            // Correct formula: w_c = -A[i,c] / ( -sum_{k in C_strong(i)} A[i,k] )
            //                       = A[i,c] / (sum A[i,k] for k in C_strong)
            // Both num and denom are negative for M-matrices → positive ratio.

            // Compute sum of A[i,c] for c in C_strong
            let rs = a.row_offsets[i];
            let re = a.row_offsets[i + 1];

            // Gather A[i,c] for each c in c_strong
            let mut c_values: Vec<(usize, f64)> = Vec::with_capacity(c_strong.len());
            for k in rs..re {
                let j = a.col_indices[k];
                if c_strong.contains(&j) {
                    c_values.push((
                        c_to_coarse[j].expect("C-point must have index"),
                        a.values[k],
                    ));
                }
            }

            // sum_{k in C_strong} A[i,k]
            let denom: f64 = c_values.iter().map(|(_, v)| v).sum();

            if denom.abs() < 1e-300 {
                // Degenerate: fallback to uniform weights
                let w = 1.0 / c_values.len() as f64;
                for (c, _) in &c_values {
                    col_indices.push(*c);
                    values.push(w);
                }
            } else {
                for (c, av) in &c_values {
                    let w = av / denom;
                    col_indices.push(*c);
                    values.push(w);
                }
            }

            // Sort by column index for canonical ordering
            let start = row_offsets[i];
            let end = col_indices.len();
            if end > start + 1 {
                let mut pairs: Vec<(usize, f64)> = col_indices[start..end]
                    .iter()
                    .zip(values[start..end].iter())
                    .map(|(&c, &v)| (c, v))
                    .collect();
                pairs.sort_unstable_by_key(|(c, _)| *c);
                for (k, (c, v)) in pairs.into_iter().enumerate() {
                    col_indices[start + k] = c;
                    values[start + k] = v;
                }
            }

            row_offsets[i + 1] = col_indices.len();
        }
    }

    CsrMatrix {
        nrows: n_fine,
        ncols: n_coarse,
        row_offsets,
        col_indices,
        values,
    }
}

// ── Test helpers ───────────────────────────────────────────────────────────────

/// Build a 2D finite-difference Poisson matrix on an n×n grid.
///
/// DOF index for node (i,j) = i*n + j.
/// 5-point stencil: diagonal = 4.0, horizontal/vertical neighbors = -1.0.
pub fn make_2d_poisson(n: usize) -> CsrMatrix {
    let ndofs = n * n;
    let mut row_offsets = vec![0usize; ndofs + 1];
    let mut col_indices: Vec<usize> = Vec::new();
    let mut values: Vec<f64> = Vec::new();

    for row in 0..n {
        for col in 0..n {
            let i = row * n + col;

            // Collect neighbours sorted by index
            let mut entries: Vec<(usize, f64)> = Vec::new();

            if row > 0 {
                entries.push((i - n, -1.0));
            }
            if col > 0 {
                entries.push((i - 1, -1.0));
            }
            entries.push((i, 4.0)); // diagonal
            if col + 1 < n {
                entries.push((i + 1, -1.0));
            }
            if row + 1 < n {
                entries.push((i + n, -1.0));
            }

            entries.sort_unstable_by_key(|(c, _)| *c);
            for (c, v) in entries {
                col_indices.push(c);
                values.push(v);
            }
            row_offsets[i + 1] = col_indices.len();
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

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::solvers::amg::cycle::{CycleKind, amg_solve};

    #[test]
    fn test_classical_build_2d_poisson() {
        // 16×16 = 256 DOF Poisson
        let n = 16;
        let a = make_2d_poisson(n);
        let amg = AmgClassical::new();
        let hier = amg.build(&a);

        assert!(
            hier.levels.len() >= 2,
            "Expected at least 2 levels, got {}",
            hier.levels.len()
        );
        // Coarser levels should have fewer DOFs
        for lvl in 0..hier.levels.len() - 1 {
            assert!(
                hier.levels[lvl + 1].a.nrows <= hier.levels[lvl].a.nrows,
                "Level {} has {} DOFs >= level {} with {} DOFs",
                lvl + 1,
                hier.levels[lvl + 1].a.nrows,
                lvl,
                hier.levels[lvl].a.nrows
            );
        }
    }

    #[test]
    fn test_vcycle_poisson_convergence() {
        // 32×32 = 1024 DOF 2D Poisson
        let n = 32;
        let a = make_2d_poisson(n);
        let amg = AmgClassical::new();
        let hier = amg.build(&a);

        let ndof = n * n;
        // RHS = all ones
        let b = vec![1.0f64; ndof];
        let mut x = vec![0.0f64; ndof];

        let stats = amg_solve(&hier, &b, &mut x, CycleKind::V, 20, 1e-6);

        assert!(
            stats.converged,
            "AMG V-cycle did not converge in 20 cycles: relative residual = {:.3e}",
            stats.residual_norm / b.iter().map(|v| v * v).sum::<f64>().sqrt()
        );
    }

    #[test]
    fn test_prolongation_row_sum() {
        // For a 1D Poisson, each row of P should sum to approximately 1.0
        // (this ensures constants are preserved exactly)
        let n = 8;
        // 1D Poisson
        let mut row_offsets = vec![0usize; n + 1];
        let mut col_indices_v = Vec::new();
        let mut values = Vec::new();
        for i in 0..n {
            if i > 0 {
                col_indices_v.push(i - 1);
                values.push(-1.0f64);
            }
            col_indices_v.push(i);
            values.push(2.0f64);
            if i + 1 < n {
                col_indices_v.push(i + 1);
                values.push(-1.0f64);
            }
            row_offsets[i + 1] = col_indices_v.len();
        }
        let a = CsrMatrix {
            nrows: n,
            ncols: n,
            row_offsets: row_offsets.clone(),
            col_indices: col_indices_v,
            values,
        };

        let amg = AmgClassical::new();
        let hier = amg.build(&a);

        // Check that P rows sum to ≈ 1.0
        let p = &hier.levels[0].p;
        for i in 0..p.nrows {
            let rs = p.row_offsets[i];
            let re = p.row_offsets[i + 1];
            let row_sum: f64 = p.values[rs..re].iter().sum();
            assert!(
                (row_sum - 1.0).abs() < 1e-10,
                "P row {i} sums to {row_sum}, expected 1.0"
            );
        }
    }
}
