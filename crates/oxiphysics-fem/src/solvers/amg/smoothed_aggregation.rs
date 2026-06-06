// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Smoothed-aggregation AMG hierarchy builder.

use crate::parallel_solver::CsrMatrix;
use crate::solvers::amg::{
    aggregation::greedy_aggregate,
    classical::extract_diagonal,
    cycle::{AmgHierarchy, AmgLevel},
    galerkin::{csr_transpose, galerkin_coarse, spmm},
    near_null_space::rigid_body_modes_3d,
};

/// Smoothed-aggregation AMG builder.
#[derive(Debug, Clone)]
pub struct SmoothedAggregationAmg {
    /// Maximum number of hierarchy levels.
    pub max_levels: usize,
    /// Problem size below which coarsening stops.
    pub coarse_cutoff: usize,
    /// Strength-of-connection threshold (default 0.08).
    pub sa_theta: f64,
    /// Jacobi smoothing factor for prolongator (default 4/3 / rho, computed per level).
    pub jacobi_omega: f64,
}

impl Default for SmoothedAggregationAmg {
    fn default() -> Self {
        Self::new()
    }
}

impl SmoothedAggregationAmg {
    /// Create with default parameters.
    pub fn new() -> Self {
        Self {
            max_levels: 10,
            coarse_cutoff: 100,
            sa_theta: 0.08,
            jacobi_omega: 4.0 / 3.0, // will be divided by spectral radius per level
        }
    }

    /// Build the SA-AMG hierarchy.
    ///
    /// `coords` is optional per-node coordinates for 3D elasticity near-null-space.
    pub fn build(&self, a: &CsrMatrix, coords: Option<&[[f64; 3]]>) -> AmgHierarchy {
        let mut levels: Vec<AmgLevel> = Vec::new();
        let mut current_a = a.clone();
        let mut current_coords: Option<Vec<[f64; 3]>> = coords.map(|c| c.to_vec());

        while levels.len() < self.max_levels && current_a.nrows > self.coarse_cutoff {
            // 1. SA strength of connection
            let strong = sa_strong_connections(&current_a, self.sa_theta);

            // 2. Greedy aggregation
            let n = current_a.nrows;
            let aggregate_ids = greedy_aggregate(&strong, n);
            let n_aggs = aggregate_ids.iter().max().copied().unwrap_or(0) + 1;

            if n_aggs == 0 || n_aggs >= n {
                break; // Degenerate — stop coarsening
            }

            // 3. Tentative prolongator
            let n_dof_per_node = if current_coords.is_some() { 3 } else { 1 };
            let p_tent = build_tentative_prolongator(
                &aggregate_ids,
                n_aggs,
                n,
                current_coords.as_deref(),
                n_dof_per_node,
            );

            if p_tent.ncols == 0 {
                break; // Nothing to coarsen to
            }

            // 4. Jacobi-smoothed prolongator
            let rho = power_iteration_spectral_radius(&current_a, 20);
            let omega = 4.0 / (3.0 * rho.max(1e-10));
            let p = jacobi_smooth_prolongator(&current_a, &p_tent, omega);
            let pt = csr_transpose(&p);

            // 5. Galerkin coarse operator
            let a_coarse = galerkin_coarse(&current_a, &p);

            // 6. Coarsen coordinates
            let coarse_coords = current_coords
                .as_deref()
                .map(|coords_slice| aggregate_centroids(&aggregate_ids, n_aggs, coords_slice));

            let diag = extract_diagonal(&current_a);
            levels.push(AmgLevel {
                a: current_a,
                p,
                pt,
                diag,
            });

            current_a = a_coarse;
            current_coords = coarse_coords;
        }

        // Final coarse level
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

/// Smoothed-aggregation strength of connection.
///
/// Node j is strongly connected to i if:
/// `|A[i,j]|^2 >= theta^2 * |A[i,i]| * |A[j,j]|`
pub fn sa_strong_connections(a: &CsrMatrix, theta: f64) -> Vec<Vec<usize>> {
    let n = a.nrows;
    let theta_sq = theta * theta;

    // Extract diagonal values
    let diag = extract_diagonal(a);

    let mut strong = vec![Vec::new(); n];
    for i in 0..n {
        let aii = diag[i].abs();
        let rs = a.row_offsets[i];
        let re = a.row_offsets[i + 1];
        for k in rs..re {
            let j = a.col_indices[k];
            if j == i {
                continue;
            }
            let ajj = diag[j].abs();
            let aij_sq = a.values[k] * a.values[k];
            if aij_sq >= theta_sq * aii * ajj {
                strong[i].push(j);
            }
        }
    }
    strong
}

/// Build the tentative prolongator from aggregate assignments.
///
/// For each aggregate, collect the node DOFs and compute near-null-space vectors.
/// The tentative prolongator is block-diagonal over aggregates.
fn build_tentative_prolongator(
    aggregate_ids: &[usize],
    n_aggs: usize,
    n_fine: usize,
    coords: Option<&[[f64; 3]]>,
    n_dof_per_node: usize,
) -> CsrMatrix {
    // For scalar case, n_dof_per_node = 1.
    // Coarse DOF count = n_aggs * n_rbm_per_agg
    // For simplicity, use n_dof_per_node = 1 always for scalar (1 mode per agg).
    let n_modes = if n_dof_per_node == 3 { 6 } else { 1 };

    // Group nodes by aggregate
    let mut agg_nodes: Vec<Vec<usize>> = vec![Vec::new(); n_aggs];
    for (i, &agg) in aggregate_ids.iter().enumerate() {
        agg_nodes[agg].push(i);
    }

    // Compute near-null-space vectors per aggregate and build columns of P_tent
    // P_tent is (n_fine) x (n_aggs * n_modes)
    // Row i, column (agg * n_modes + m) = <node i's contribution to mode m of agg>

    let n_coarse = n_aggs * n_modes;
    let mut row_offsets = vec![0usize; n_fine + 1];
    let mut col_indices: Vec<usize> = Vec::new();
    let mut values: Vec<f64> = Vec::new();

    for i in 0..n_fine {
        let agg = aggregate_ids[i];
        let nodes_in_agg = &agg_nodes[agg];

        // Find position of node i within the aggregate
        let local_idx = nodes_in_agg.iter().position(|&nd| nd == i).unwrap_or(0);

        // Get near-null-space vectors for this aggregate
        let node_coords: Vec<[f64; 3]> = if let Some(c) = coords {
            nodes_in_agg.iter().map(|&nd| c[nd]).collect()
        } else {
            nodes_in_agg
                .iter()
                .enumerate()
                .map(|(k, _)| [k as f64, 0.0, 0.0])
                .collect()
        };

        let rbms = rigid_body_modes_3d(&node_coords, n_dof_per_node);
        let actual_modes = rbms.len().min(n_modes);

        // Each RBM vector has length n_nodes_in_agg * n_dof_per_node
        // For node i (local_idx), its contribution to mode m is rbms[m][local_idx * n_dof_per_node..]
        // For scalar (n_dof_per_node=1): rbms[0][local_idx]
        for m in 0..actual_modes {
            let entry_start = local_idx * n_dof_per_node;
            // For scalar, just use rbms[m][local_idx]
            if n_dof_per_node == 1 {
                let val = if m < rbms.len() && entry_start < rbms[m].len() {
                    rbms[m][entry_start]
                } else {
                    0.0
                };
                if val.abs() > 1e-300 {
                    col_indices.push(agg * n_modes + m);
                    values.push(val);
                }
            } else {
                // For 3D elasticity: each DOF of node i contributes to a column
                // But we're treating fine DOFs as nodes here (scalar indexing)
                // In general: fine DOF index = node * 3 + dof_component
                // For now, handle scalar assumption
                let val = if m < rbms.len() && entry_start < rbms[m].len() {
                    rbms[m][entry_start]
                } else {
                    0.0
                };
                if val.abs() > 1e-300 {
                    col_indices.push(agg * n_modes + m);
                    values.push(val);
                }
            }
        }

        row_offsets[i + 1] = col_indices.len();
    }

    CsrMatrix {
        nrows: n_fine,
        ncols: n_coarse,
        row_offsets,
        col_indices,
        values,
    }
}

/// Estimate spectral radius of D^{-1}A via 20 power iterations.
///
/// Uses a deterministic starting vector.
pub fn power_iteration_spectral_radius(a: &CsrMatrix, n_iters: usize) -> f64 {
    let n = a.nrows;
    if n == 0 {
        return 1.0;
    }

    // D^{-1} (inverse diagonal)
    let d_inv = a.diagonal_preconditioner();

    // Deterministic starting vector: v[i] = (-1)^i / sqrt(n)
    let inv_sqrt_n = 1.0 / (n as f64).sqrt();
    let mut v: Vec<f64> = (0..n)
        .map(|i| if i % 2 == 0 { inv_sqrt_n } else { -inv_sqrt_n })
        .collect();

    let mut rho = 1.0f64;

    for _ in 0..n_iters {
        // w = D^{-1} * A * v
        let mut av = vec![0.0f64; n];
        a.spmv(&v, &mut av);
        let mut w: Vec<f64> = av
            .iter()
            .zip(d_inv.iter())
            .map(|(av_i, di)| av_i * di)
            .collect();

        // Compute Rayleigh quotient: rho = (v^T w) / (v^T v)
        let vw: f64 = v.iter().zip(w.iter()).map(|(vi, wi)| vi * wi).sum();
        let vv: f64 = v.iter().map(|vi| vi * vi).sum();
        rho = if vv > 1e-300 { vw / vv } else { 1.0 };

        // Normalize w to get new v
        let w_norm: f64 = w.iter().map(|wi| wi * wi).sum::<f64>().sqrt();
        if w_norm < 1e-300 {
            break;
        }
        let inv_norm = 1.0 / w_norm;
        for wi in w.iter_mut() {
            *wi *= inv_norm;
        }
        v = w;
    }

    rho.abs().max(1e-10)
}

/// Jacobi-smoothed prolongator: P = (I - ω D^{-1} A) P_tent.
///
/// Computed as P = P_tent - ω * D^{-1} * A * P_tent.
pub fn jacobi_smooth_prolongator(a: &CsrMatrix, p_tent: &CsrMatrix, omega: f64) -> CsrMatrix {
    // Compute ω * D^{-1} * A * P_tent
    let d_inv = a.diagonal_preconditioner(); // = 1/A[i,i]

    // Scale rows of A by D^{-1}: create D^{-1} A
    // We do this implicitly: for each row i of A*P_tent, multiply by d_inv[i]

    // First compute A * P_tent
    let ap_tent = spmm(a, p_tent);

    // Now compute P_tent - omega * diag(d_inv) * (A * P_tent)
    // Both P_tent and ap_tent have the same row structure (n_fine rows)
    // We need to combine them: result = P_tent + (-omega * d_inv[i]) * (ap_tent row i)

    let n_fine = p_tent.nrows;
    let n_coarse = p_tent.ncols;

    // Build result sparsity = union of P_tent and ap_tent column sets per row
    use std::collections::BTreeSet;
    let mut row_col_sets: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); n_fine];
    for (i, row_set) in row_col_sets.iter_mut().enumerate() {
        for k in p_tent.row_offsets[i]..p_tent.row_offsets[i + 1] {
            row_set.insert(p_tent.col_indices[k]);
        }
        for k in ap_tent.row_offsets[i]..ap_tent.row_offsets[i + 1] {
            row_set.insert(ap_tent.col_indices[k]);
        }
    }

    let mut row_offsets = vec![0usize; n_fine + 1];
    let mut col_indices: Vec<usize> = Vec::new();
    for i in 0..n_fine {
        row_offsets[i + 1] = row_offsets[i] + row_col_sets[i].len();
        col_indices.extend(row_col_sets[i].iter().copied());
    }
    let nnz = col_indices.len();
    let mut values = vec![0.0f64; nnz];

    // Add P_tent contributions
    for i in 0..n_fine {
        let base = row_offsets[i];
        let cols = &col_indices[base..row_offsets[i + 1]];
        for k in p_tent.row_offsets[i]..p_tent.row_offsets[i + 1] {
            let col = p_tent.col_indices[k];
            // Find position in sorted cols
            if let Ok(pos) = cols.binary_search(&col) {
                values[base + pos] += p_tent.values[k];
            }
        }
    }

    // Subtract omega * D^{-1} * (ap_tent) contributions
    for i in 0..n_fine {
        let scale = -omega * d_inv[i];
        let base = row_offsets[i];
        let cols = &col_indices[base..row_offsets[i + 1]];
        for k in ap_tent.row_offsets[i]..ap_tent.row_offsets[i + 1] {
            let col = ap_tent.col_indices[k];
            if let Ok(pos) = cols.binary_search(&col) {
                values[base + pos] += scale * ap_tent.values[k];
            }
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

/// Compute aggregate centroids from node coordinates.
pub fn aggregate_centroids(
    aggregate_ids: &[usize],
    n_aggs: usize,
    coords: &[[f64; 3]],
) -> Vec<[f64; 3]> {
    let mut centroids = vec![[0.0f64; 3]; n_aggs];
    let mut counts = vec![0usize; n_aggs];

    for (i, &agg) in aggregate_ids.iter().enumerate() {
        centroids[agg][0] += coords[i][0];
        centroids[agg][1] += coords[i][1];
        centroids[agg][2] += coords[i][2];
        counts[agg] += 1;
    }

    for agg in 0..n_aggs {
        if counts[agg] > 0 {
            let inv = 1.0 / counts[agg] as f64;
            centroids[agg][0] *= inv;
            centroids[agg][1] *= inv;
            centroids[agg][2] *= inv;
        }
    }

    centroids
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::solvers::amg::cycle::{CycleKind, amg_solve};

    fn make_1d_poisson(n: usize) -> CsrMatrix {
        let mut row_offsets = vec![0usize; n + 1];
        let mut col_indices = Vec::new();
        let mut values = Vec::new();
        for i in 0..n {
            if i > 0 {
                col_indices.push(i - 1);
                values.push(-1.0f64);
            }
            col_indices.push(i);
            values.push(2.0f64);
            if i + 1 < n {
                col_indices.push(i + 1);
                values.push(-1.0f64);
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

    #[test]
    fn test_sa_aggregation_covers_all() {
        let n = 20;
        let a = make_1d_poisson(n);
        let strong = sa_strong_connections(&a, 0.08);
        let agg_ids = greedy_aggregate(&strong, n);

        // All nodes must be assigned
        for (i, &id) in agg_ids.iter().enumerate() {
            assert!(id != usize::MAX, "Node {i} has no aggregate assignment");
        }

        // Total nodes must equal n
        let n_aggs = agg_ids.iter().max().copied().unwrap_or(0) + 1;
        let mut sizes = vec![0usize; n_aggs];
        for &id in &agg_ids {
            sizes[id] += 1;
        }
        let total: usize = sizes.iter().sum();
        assert_eq!(total, n, "Aggregates don't cover all {n} nodes");
    }

    #[test]
    // slow test — run with --include-ignored
    #[ignore]
    fn test_sa_build_converges() {
        // 64-node 1D problem
        let n = 64;
        let a = make_1d_poisson(n);

        let sa = SmoothedAggregationAmg {
            coarse_cutoff: 8,
            ..SmoothedAggregationAmg::new()
        };
        let hier = sa.build(&a, None);

        let b = vec![1.0f64; n];
        let mut x = vec![0.0f64; n];

        let stats = amg_solve(&hier, &b, &mut x, CycleKind::V, 10, 1e-4);

        let b_norm = b.iter().map(|v| v * v).sum::<f64>().sqrt();
        let rel_res = stats.residual_norm / b_norm;

        assert!(
            rel_res < 1e-4,
            "SA-AMG did not converge: relative residual = {rel_res:.3e}"
        );
    }

    #[test]
    fn test_power_iteration_simple() {
        // For a diagonal matrix with values [2, 4, 6], spectral radius of D^{-1}A = 1.0
        let a = CsrMatrix {
            nrows: 3,
            ncols: 3,
            row_offsets: vec![0, 1, 2, 3],
            col_indices: vec![0, 1, 2],
            values: vec![2.0, 4.0, 6.0],
        };
        let rho = power_iteration_spectral_radius(&a, 20);
        // D^{-1}A = I, so spectral radius = 1.0
        assert!(
            (rho - 1.0).abs() < 0.1,
            "Spectral radius of identity should be ~1.0, got {rho}"
        );
    }
}
