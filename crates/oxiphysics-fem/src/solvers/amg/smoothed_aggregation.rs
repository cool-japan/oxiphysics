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
///
/// Two kinds of problems are supported:
///
/// * **Scalar** operators (one DOF per node; Poisson-type): the near-null space
///   is the constant vector.
/// * **3-D elasticity** (three displacement DOFs per node, ordered
///   `(u_x, u_y, u_z)` per node): the near-null space is the six rigid-body
///   modes (three translations, three rotations).  Selected by passing the node
///   coordinates to [`build`](Self::build).
///
/// The hierarchy follows Vaněk, Mandel and Brezina (1996).  At every level the
/// nodes are aggregated on the node graph; for each aggregate the restriction
/// of the near-null-space matrix `B` to the aggregate's DOFs is factored
/// `B_agg = Q R` (thin QR with rank detection).  `Q` becomes the aggregate's
/// block of the tentative prolongator, and the factor `R` is passed down as the
/// near-null space of the coarse level, so that `P_tent B_coarse = B` holds
/// exactly and every level reproduces the rigid-body modes.  Aggregates may
/// contribute fewer than six coarse DOFs (a single node carries only the three
/// translations, two nodes five modes), so the coarse levels have a variable
/// number of DOFs per node.
#[derive(Debug, Clone)]
pub struct SmoothedAggregationAmg {
    /// Maximum number of hierarchy levels.
    pub max_levels: usize,
    /// Problem size below which coarsening stops.
    pub coarse_cutoff: usize,
    /// Strength-of-connection threshold (default 0.08).
    pub sa_theta: f64,
    /// Numerator of the Jacobi prolongator-smoothing factor (default 4/3); the
    /// damping used at each level is `jacobi_omega / ρ(D⁻¹A)` with the spectral
    /// radius estimated per level.
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
    /// `coords` selects the near-null space.  With `None`, or if
    /// `3 * coords.len() != a.nrows`, the operator is treated as a scalar one
    /// (constant near-null space, one DOF per node).  With one coordinate per
    /// node and `a.nrows == 3 * coords.len()` the operator is a 3-D elasticity
    /// matrix with DOFs ordered `(u_x, u_y, u_z)` per node, and the six
    /// rigid-body modes generated from the coordinates are the near-null space.
    ///
    /// Nodes that have no strong connection (for example DOFs eliminated by
    /// Dirichlet conditions with an identity row) are not aggregated: their
    /// rows of the prolongator are zero.
    pub fn build(&self, a: &CsrMatrix, coords: Option<&[[f64; 3]]>) -> AmgHierarchy {
        let mut levels: Vec<AmgLevel> = Vec::new();
        let mut current_a = a.clone();
        let (mut node_offsets, mut null_space) = initial_node_structure(a.nrows, coords);

        while levels.len() < self.max_levels && current_a.nrows > self.coarse_cutoff {
            let n = current_a.nrows;
            let n_nodes = node_offsets.len() - 1;

            // 1. Node-level strength of connection (Frobenius norm of blocks)
            let dof_to_node = dof_to_node_map(&node_offsets);
            let strong =
                node_strong_connections(&current_a, &node_offsets, &dof_to_node, self.sa_theta);

            // 2. Greedy aggregation of the nodes that have strong connections
            let node_aggregate = aggregate_connected_nodes(&strong, n_nodes);
            let n_aggs = node_aggregate
                .iter()
                .filter(|&&agg| agg != NO_AGGREGATE)
                .max()
                .map_or(0, |&top| top + 1);
            if n_aggs == 0 {
                break; // Nothing to aggregate
            }

            // 3. Tentative prolongator and coarse near-null space (thin QR)
            let tent =
                build_tentative_prolongator(&node_aggregate, n_aggs, &node_offsets, &null_space);
            let n_coarse = tent.p.ncols;
            if n_coarse == 0 || n_coarse >= n {
                break; // Degenerate: coarsening would not reduce the problem
            }

            // 4. Jacobi-smoothed prolongator
            let rho = power_iteration_spectral_radius(&current_a, 20);
            let omega = self.jacobi_omega / rho.max(1e-10);
            let p = jacobi_smooth_prolongator(&current_a, &tent.p, omega);
            let pt = csr_transpose(&p);

            // 5. Galerkin coarse operator
            let a_coarse = galerkin_coarse(&current_a, &p);

            let diag = extract_diagonal(&current_a);
            levels.push(AmgLevel {
                a: current_a,
                p,
                pt,
                diag,
            });

            current_a = a_coarse;
            node_offsets = tent.coarse_node_offsets;
            null_space = tent.coarse_null_space;
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

/// Marker for a node that belongs to no aggregate.
const NO_AGGREGATE: usize = usize::MAX;

/// Near-null-space matrix `B` (rows = DOFs, `k` columns) stored row-major.
#[derive(Debug, Clone)]
struct NullSpace {
    /// Number of near-null-space vectors.
    k: usize,
    /// `values[dof * k + c]` is the entry of vector `c` at DOF `dof`.
    values: Vec<f64>,
}

/// Tentative prolongator together with the coarse-level node structure.
struct TentativeProlongator {
    /// `P_tent` (`n_fine × n_coarse`) with orthonormal columns per aggregate.
    p: CsrMatrix,
    /// DOF offsets of the coarse nodes (one per non-empty aggregate).
    coarse_node_offsets: Vec<usize>,
    /// Near-null space of the coarse level: the `R` factors of the QR steps.
    coarse_null_space: NullSpace,
}

/// DOF offsets of the nodes and the near-null space of the finest level.
fn initial_node_structure(n: usize, coords: Option<&[[f64; 3]]>) -> (Vec<usize>, NullSpace) {
    if let Some(c) = coords
        && !c.is_empty()
        && 3 * c.len() == n
    {
        // Rigid-body modes of the whole mesh, orthonormalised; dependent
        // modes (planar or collinear node sets) have been removed.
        let modes = rigid_body_modes_3d(c, 3);
        if !modes.is_empty() {
            let k = modes.len();
            let mut values = vec![0.0_f64; n * k];
            for (col, mode) in modes.iter().enumerate() {
                for (dof, &entry) in mode.iter().enumerate() {
                    values[dof * k + col] = entry;
                }
            }
            let node_offsets = (0..=c.len()).map(|node| 3 * node).collect();
            return (node_offsets, NullSpace { k, values });
        }
    }
    let node_offsets = (0..=n).collect();
    (
        node_offsets,
        NullSpace {
            k: 1,
            values: vec![1.0; n],
        },
    )
}

/// Node index of every DOF, from the node DOF offsets.
fn dof_to_node_map(node_offsets: &[usize]) -> Vec<usize> {
    let n_dofs = node_offsets.last().copied().unwrap_or(0);
    let mut map = vec![0usize; n_dofs];
    for (node, window) in node_offsets.windows(2).enumerate() {
        for entry in &mut map[window[0]..window[1]] {
            *entry = node;
        }
    }
    map
}

/// Node-level strength of connection for a matrix whose DOFs are grouped into
/// nodes (`node_offsets`).
///
/// Node `J` is strongly connected to node `I` if
/// `‖A_IJ‖_F² ≥ θ² ‖A_II‖_F ‖A_JJ‖_F`, the block analogue of the scalar test of
/// [`sa_strong_connections`].  Neighbour lists are sorted by node index.
fn node_strong_connections(
    a: &CsrMatrix,
    node_offsets: &[usize],
    dof_to_node: &[usize],
    theta: f64,
) -> Vec<Vec<usize>> {
    let n_nodes = node_offsets.len() - 1;
    let theta_sq = theta * theta;
    let row_range = |row: usize| a.row_offsets[row]..a.row_offsets[row + 1];

    let mut diag_block_norm = vec![0.0_f64; n_nodes];
    for (node, window) in node_offsets.windows(2).enumerate() {
        let mut sum_sq = 0.0;
        for row in window[0]..window[1] {
            for k in row_range(row) {
                if dof_to_node[a.col_indices[k]] == node {
                    sum_sq += a.values[k] * a.values[k];
                }
            }
        }
        diag_block_norm[node] = sum_sq.sqrt();
    }

    let mut strong = vec![Vec::new(); n_nodes];
    let mut block_sq = vec![0.0_f64; n_nodes];
    let mut stamp = vec![usize::MAX; n_nodes];
    let mut touched: Vec<usize> = Vec::new();
    for (node, window) in node_offsets.windows(2).enumerate() {
        for row in window[0]..window[1] {
            for k in row_range(row) {
                let other = dof_to_node[a.col_indices[k]];
                if other == node {
                    continue;
                }
                if stamp[other] != node {
                    stamp[other] = node;
                    block_sq[other] = 0.0;
                    touched.push(other);
                }
                block_sq[other] += a.values[k] * a.values[k];
            }
        }
        touched.sort_unstable();
        for &other in &touched {
            if block_sq[other] > 0.0
                && block_sq[other] >= theta_sq * diag_block_norm[node] * diag_block_norm[other]
            {
                strong[node].push(other);
            }
        }
        touched.clear();
    }
    strong
}

/// Greedy aggregation of the nodes that have at least one strong connection.
///
/// Returns the aggregate id of every node, or [`NO_AGGREGATE`] for isolated
/// nodes (no strong neighbour, typically DOFs eliminated by an identity row).
fn aggregate_connected_nodes(strong: &[Vec<usize>], n_nodes: usize) -> Vec<usize> {
    let mut compact = vec![NO_AGGREGATE; n_nodes];
    let mut active: Vec<usize> = Vec::new();
    for (node, neighbours) in strong.iter().enumerate() {
        if !neighbours.is_empty() {
            compact[node] = active.len();
            active.push(node);
        }
    }
    if active.is_empty() {
        return vec![NO_AGGREGATE; n_nodes];
    }
    let compact_strong: Vec<Vec<usize>> = active
        .iter()
        .map(|&node| {
            strong[node]
                .iter()
                .map(|&other| compact[other])
                .filter(|&c| c != NO_AGGREGATE)
                .collect()
        })
        .collect();
    let compact_aggregate = greedy_aggregate(&compact_strong, active.len());
    let mut aggregate = vec![NO_AGGREGATE; n_nodes];
    for (c, &node) in active.iter().enumerate() {
        aggregate[node] = compact_aggregate[c];
    }
    aggregate
}

/// Thin QR factorisation `B = Q R` with rank detection, by modified
/// Gram-Schmidt with one re-orthogonalisation pass.
///
/// `columns` are the `k` columns of `B` (each of length `m`).  Returns the `r`
/// orthonormal columns of `Q` (each of length `m`) and the `r × k` factor `R`
/// (`R[i][c]` is the coefficient of `Q[i]` in column `c`).  A column whose
/// component outside the span of the previous ones is below `1e-10` of its own
/// norm is dependent and contributes no `Q` column (so `r ≤ min(m, k)`).
fn thin_qr(columns: &[Vec<f64>]) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    const RANK_TOLERANCE: f64 = 1.0e-10;
    let k = columns.len();
    let mut q: Vec<Vec<f64>> = Vec::with_capacity(k);
    let mut r: Vec<Vec<f64>> = Vec::with_capacity(k);
    for (c, column) in columns.iter().enumerate() {
        let original_norm = column.iter().map(|v| v * v).sum::<f64>().sqrt();
        let mut v = column.clone();
        let mut coefficients = vec![0.0_f64; q.len()];
        for _pass in 0..2 {
            for (qi, coefficient) in q.iter().zip(coefficients.iter_mut()) {
                let projection: f64 = qi.iter().zip(v.iter()).map(|(a, b)| a * b).sum();
                for (vj, qj) in v.iter_mut().zip(qi.iter()) {
                    *vj -= projection * qj;
                }
                *coefficient += projection;
            }
        }
        for (row, coefficient) in r.iter_mut().zip(coefficients.iter()) {
            row[c] = *coefficient;
        }
        let residual = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        if original_norm > 0.0 && residual > RANK_TOLERANCE * original_norm {
            for entry in v.iter_mut() {
                *entry /= residual;
            }
            let mut row = vec![0.0_f64; k];
            row[c] = residual;
            q.push(v);
            r.push(row);
        }
    }
    (q, r)
}

/// Build the tentative prolongator from the node aggregates.
///
/// For every aggregate the rows of the near-null-space matrix `B` belonging to
/// its DOFs are factored `B_agg = Q R` ([`thin_qr`]); the orthonormal `Q` forms
/// the aggregate's block of `P_tent`, whose columns are numbered consecutively
/// over the aggregates, and the rows of `R` become the rows of the coarse
/// near-null space, so that `P_tent B_coarse = B` on the aggregated nodes.
/// Each aggregate with `r > 0` becomes one coarse node of `r` DOFs.  DOFs of
/// nodes without an aggregate have empty rows.
fn build_tentative_prolongator(
    node_aggregate: &[usize],
    n_aggs: usize,
    node_offsets: &[usize],
    null_space: &NullSpace,
) -> TentativeProlongator {
    let n_fine = node_offsets.last().copied().unwrap_or(0);
    let k = null_space.k;

    let mut agg_nodes: Vec<Vec<usize>> = vec![Vec::new(); n_aggs];
    for (node, &agg) in node_aggregate.iter().enumerate() {
        if agg != NO_AGGREGATE {
            agg_nodes[agg].push(node);
        }
    }

    let mut agg_dofs: Vec<Vec<usize>> = Vec::with_capacity(n_aggs);
    let mut agg_q: Vec<Vec<Vec<f64>>> = Vec::with_capacity(n_aggs);
    let mut agg_first_col: Vec<usize> = Vec::with_capacity(n_aggs);
    let mut coarse_node_offsets: Vec<usize> = vec![0];
    let mut coarse_values: Vec<f64> = Vec::new();
    let mut n_coarse = 0usize;
    for nodes in &agg_nodes {
        let dofs: Vec<usize> = nodes
            .iter()
            .flat_map(|&node| node_offsets[node]..node_offsets[node + 1])
            .collect();
        let columns: Vec<Vec<f64>> = (0..k)
            .map(|c| {
                dofs.iter()
                    .map(|&dof| null_space.values[dof * k + c])
                    .collect()
            })
            .collect();
        let (q, r) = thin_qr(&columns);
        agg_first_col.push(n_coarse);
        if !q.is_empty() {
            n_coarse += q.len();
            coarse_node_offsets.push(n_coarse);
            for row in &r {
                coarse_values.extend_from_slice(row);
            }
        }
        agg_dofs.push(dofs);
        agg_q.push(q);
    }

    let mut row_offsets = vec![0usize; n_fine + 1];
    for (dofs, q) in agg_dofs.iter().zip(agg_q.iter()) {
        for &dof in dofs {
            row_offsets[dof + 1] = q.len();
        }
    }
    for i in 0..n_fine {
        row_offsets[i + 1] += row_offsets[i];
    }
    let nnz = row_offsets[n_fine];
    let mut col_indices = vec![0usize; nnz];
    let mut values = vec![0.0_f64; nnz];
    for ((dofs, q), &first_col) in agg_dofs.iter().zip(agg_q.iter()).zip(agg_first_col.iter()) {
        for (local, &dof) in dofs.iter().enumerate() {
            let base = row_offsets[dof];
            for (j, q_col) in q.iter().enumerate() {
                col_indices[base + j] = first_col + j;
                values[base + j] = q_col[local];
            }
        }
    }

    TentativeProlongator {
        p: CsrMatrix {
            nrows: n_fine,
            ncols: n_coarse,
            row_offsets,
            col_indices,
            values,
        },
        coarse_node_offsets,
        coarse_null_space: NullSpace {
            k,
            values: coarse_values,
        },
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
///
/// `aggregate_ids[i]` is the aggregate of *node* `i` and `coords[i]` its
/// position (one entry per node, not per DOF).  Nodes with an aggregate id
/// `>= n_aggs` (for example unaggregated nodes) are ignored, and so are
/// surplus entries when the two slices differ in length.
pub fn aggregate_centroids(
    aggregate_ids: &[usize],
    n_aggs: usize,
    coords: &[[f64; 3]],
) -> Vec<[f64; 3]> {
    let mut centroids = vec![[0.0f64; 3]; n_aggs];
    let mut counts = vec![0usize; n_aggs];

    for (&agg, coord) in aggregate_ids.iter().zip(coords.iter()) {
        if agg >= n_aggs {
            continue;
        }
        centroids[agg][0] += coord[0];
        centroids[agg][1] += coord[1];
        centroids[agg][2] += coord[2];
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
    fn test_aggregate_centroids_ignores_unaggregated_and_surplus_nodes() {
        let coords = [
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 3.0, 0.0],
            [9.0, 9.0, 9.0],
        ];
        // Node 3 is unaggregated; the id list is longer than the coordinates.
        let ids = [0, 0, 1, NO_AGGREGATE, 1, 0];
        let centroids = aggregate_centroids(&ids, 2, &coords);
        assert_eq!(centroids[0], [1.0, 0.0, 0.0]);
        assert_eq!(centroids[1], [0.0, 3.0, 0.0]);
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

    // ── Block (elasticity) near-null-space tests ─────────────────────────────

    use crate::assembly::{apply_dirichlet_bc, assemble_stiffness};
    use crate::constitutive::LinearElasticMaterial;
    use crate::mesh::TetrahedralMesh;
    use crate::parallel_solver::{ParallelPcgSolver, PcgWithPrecond};
    use crate::solvers::CholeskyDense;
    use crate::solvers::amg::preconditioner::AmgPreconditioner;

    /// Convert the assembler's CSR format to the AMG module's CSR format.
    fn to_parallel_csr(k: &crate::sparse::CsrMatrix) -> CsrMatrix {
        CsrMatrix {
            nrows: k.nrows,
            ncols: k.ncols,
            row_offsets: k.row_ptr.clone(),
            col_indices: k.col_indices.clone(),
            values: k.values.clone(),
        }
    }

    /// Linear-elasticity stiffness matrix of the unit cube meshed with
    /// `cells³` hexahedral cells (5 tetrahedra each), with node coordinates.
    /// If `clamp_bottom`, all DOFs of the `z = 0` face are eliminated by
    /// identity rows (the matrix keeps its size `3 · n_nodes`).
    fn cube_elasticity(cells: usize, clamp_bottom: bool) -> (CsrMatrix, Vec<[f64; 3]>, Vec<usize>) {
        let mesh = TetrahedralMesh::generate_beam(1.0, 1.0, 1.0, cells, cells, cells);
        let material = LinearElasticMaterial::new(1.0e3, 0.3);
        let mut k = assemble_stiffness(&mesh, &material);
        let coords: Vec<[f64; 3]> = mesh.nodes.iter().map(|p| [p.x, p.y, p.z]).collect();
        let mut fixed: Vec<usize> = Vec::new();
        if clamp_bottom {
            for (node, c) in coords.iter().enumerate() {
                if c[2] == 0.0 {
                    fixed.extend([3 * node, 3 * node + 1, 3 * node + 2]);
                }
            }
            let mut f = vec![0.0_f64; k.nrows];
            apply_dirichlet_bc(&mut k, &mut f, &fixed).expect("fixed DOFs are in range");
        }
        (to_parallel_csr(&k), coords, fixed)
    }

    fn dense_row(a: &CsrMatrix, i: usize) -> Vec<(usize, f64)> {
        (a.row_offsets[i]..a.row_offsets[i + 1])
            .map(|k| (a.col_indices[k], a.values[k]))
            .collect()
    }

    /// Least-squares residual `‖b − P c‖ / ‖b‖` of the best `c` (normal equations).
    fn range_residual(p: &CsrMatrix, b: &[f64]) -> f64 {
        let nc = p.ncols;
        let mut gram = vec![0.0_f64; nc * nc];
        let mut rhs = vec![0.0_f64; nc];
        for (i, &b_i) in b.iter().enumerate().take(p.nrows) {
            let row = dense_row(p, i);
            for &(c1, v1) in &row {
                rhs[c1] += v1 * b_i;
                for &(c2, v2) in &row {
                    gram[c1 * nc + c2] += v1 * v2;
                }
            }
        }
        let chol = CholeskyDense::new(nc);
        let factor = chol.factorize(&gram).expect("P has full column rank");
        let coeff = chol.solve(&factor, &rhs);
        let mut res_sq = 0.0;
        for (i, &b_i) in b.iter().enumerate().take(p.nrows) {
            let fit: f64 = dense_row(p, i).iter().map(|&(c, v)| v * coeff[c]).sum();
            res_sq += (b_i - fit).powi(2);
        }
        res_sq.sqrt() / b.iter().map(|v| v * v).sum::<f64>().sqrt()
    }

    #[test]
    fn test_thin_qr_reproduces_columns_and_drops_dependent_ones() {
        // Columns: a, b, a + 2b (dependent), a zero column, and a new direction.
        let a = vec![1.0, 2.0, 0.0, 1.0];
        let b = vec![0.0, 1.0, 1.0, -1.0];
        let dep: Vec<f64> = a.iter().zip(b.iter()).map(|(x, y)| x + 2.0 * y).collect();
        let zero = vec![0.0; 4];
        let fresh = vec![1.0, 0.0, 0.0, 0.0];
        let columns = vec![a, b, dep, zero, fresh];
        let (q, r) = thin_qr(&columns);
        assert_eq!(q.len(), 3, "rank of the column set is 3");
        assert_eq!(r.len(), 3);
        for i in 0..3 {
            for j in 0..3 {
                let d: f64 = q[i].iter().zip(q[j].iter()).map(|(x, y)| x * y).sum();
                let want = if i == j { 1.0 } else { 0.0 };
                assert!((d - want).abs() < 1e-14, "Q^T Q[{i}][{j}] = {d}");
            }
        }
        for (c, column) in columns.iter().enumerate() {
            for row in 0..4 {
                let rebuilt: f64 = (0..3).map(|i| q[i][row] * r[i][c]).sum();
                assert!(
                    (rebuilt - column[row]).abs() < 1e-13,
                    "column {c}, row {row}: {rebuilt} vs {}",
                    column[row]
                );
            }
        }
        // R is upper triangular: q_i does not appear in columns before it.
        assert_eq!(r[2][0], 0.0);
        assert_eq!(r[2][1], 0.0);
        assert_eq!(r[1][0], 0.0);
    }

    #[test]
    fn test_node_strong_connections_use_block_frobenius_norms() {
        // Three nodes with 2 DOFs each: nodes 0-1 coupled strongly, 1-2 weakly.
        let mut triplets: Vec<(usize, usize, f64)> = Vec::new();
        for dof in 0..6 {
            triplets.push((dof, dof, 4.0));
        }
        for (r, c, v) in [(0, 2, -1.0), (1, 3, -1.0), (2, 0, -1.0), (3, 1, -1.0)] {
            triplets.push((r, c, v));
        }
        for (r, c, v) in [(2, 4, -0.01), (3, 5, -0.01), (4, 2, -0.01), (5, 3, -0.01)] {
            triplets.push((r, c, v));
        }
        let mut row_offsets = vec![0usize; 7];
        let mut col_indices = Vec::new();
        let mut values = Vec::new();
        for row in 0..6 {
            let mut entries: Vec<(usize, f64)> = triplets
                .iter()
                .filter(|t| t.0 == row)
                .map(|t| (t.1, t.2))
                .collect();
            entries.sort_by_key(|e| e.0);
            for (c, v) in entries {
                col_indices.push(c);
                values.push(v);
            }
            row_offsets[row + 1] = col_indices.len();
        }
        let a = CsrMatrix {
            nrows: 6,
            ncols: 6,
            row_offsets,
            col_indices,
            values,
        };
        let node_offsets = vec![0, 2, 4, 6];
        let dof_to_node = dof_to_node_map(&node_offsets);
        assert_eq!(dof_to_node, vec![0, 0, 1, 1, 2, 2]);
        let strong = node_strong_connections(&a, &node_offsets, &dof_to_node, 0.08);
        // ||A_01||_F^2 = 2, ||A_00||_F = ||A_11||_F = ||A_22||_F = sqrt(32):
        // 2 >= 0.0064 * sqrt(32) * sqrt(32); ||A_12||_F^2 = 2e-4 < 0.0064 * 32.
        assert_eq!(strong, vec![vec![1], vec![0], vec![]]);
        let aggregate = aggregate_connected_nodes(&strong, 3);
        assert_eq!(aggregate[0], aggregate[1]);
        assert_eq!(
            aggregate[2], NO_AGGREGATE,
            "an isolated node is not aggregated"
        );
    }

    #[test]
    fn test_tentative_prolongator_reproduces_the_near_null_space() {
        // 5 nodes with 3 DOFs, rigid-body modes from coordinates, two aggregates.
        let coords: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.2],
            [0.3, 0.2, 1.0],
            [1.2, 1.1, 0.9],
        ];
        let (node_offsets, null_space) = initial_node_structure(15, Some(&coords));
        assert_eq!(null_space.k, 6);
        let node_aggregate = vec![0, 0, 0, 1, 1];
        let tent = build_tentative_prolongator(&node_aggregate, 2, &node_offsets, &null_space);
        // Aggregate 0 has three non-collinear nodes (rank 6), aggregate 1 two
        // nodes (rank 5: the rotation about their connecting line vanishes).
        assert_eq!(tent.p.nrows, 15);
        assert_eq!(tent.p.ncols, 11);
        assert_eq!(tent.coarse_node_offsets, vec![0, 6, 11]);
        assert_eq!(tent.coarse_null_space.k, 6);
        assert_eq!(tent.coarse_null_space.values.len(), 11 * 6);
        // P_tent B_coarse = B exactly.
        for dof in 0..15 {
            for c in 0..6 {
                let mut sum = 0.0;
                for (col, value) in dense_row(&tent.p, dof) {
                    sum += value * tent.coarse_null_space.values[col * 6 + c];
                }
                assert!(
                    (sum - null_space.values[dof * 6 + c]).abs() < 1e-13,
                    "dof {dof}, mode {c}: {sum} vs {}",
                    null_space.values[dof * 6 + c]
                );
            }
        }
        // Orthonormal columns: P_tent^T P_tent = I.
        for c1 in 0..11 {
            for c2 in 0..11 {
                let mut dot = 0.0;
                for dof in 0..15 {
                    let row = dense_row(&tent.p, dof);
                    let v1 = row.iter().find(|e| e.0 == c1).map_or(0.0, |e| e.1);
                    let v2 = row.iter().find(|e| e.0 == c2).map_or(0.0, |e| e.1);
                    dot += v1 * v2;
                }
                let want = if c1 == c2 { 1.0 } else { 0.0 };
                assert!((dot - want).abs() < 1e-13, "PtP[{c1}][{c2}] = {dot}");
            }
        }
    }

    #[test]
    fn test_scalar_tentative_prolongator_is_normalised_piecewise_constant() {
        let (node_offsets, null_space) = initial_node_structure(7, None);
        assert_eq!(null_space.k, 1);
        let node_aggregate = vec![0, 0, 0, 1, 1, NO_AGGREGATE, 2];
        let tent = build_tentative_prolongator(&node_aggregate, 3, &node_offsets, &null_space);
        assert_eq!(tent.p.ncols, 3);
        assert_eq!(tent.coarse_node_offsets, vec![0, 1, 2, 3]);
        let s3 = 1.0 / 3.0_f64.sqrt();
        let s2 = 1.0 / 2.0_f64.sqrt();
        let expected = [
            vec![(0, s3)],
            vec![(0, s3)],
            vec![(0, s3)],
            vec![(1, s2)],
            vec![(1, s2)],
            vec![],
            vec![(2, 1.0)],
        ];
        for (i, want) in expected.iter().enumerate() {
            let got = dense_row(&tent.p, i);
            assert_eq!(got.len(), want.len(), "row {i}");
            for (g, w) in got.iter().zip(want.iter()) {
                assert_eq!(g.0, w.0);
                assert!((g.1 - w.1).abs() < 1e-14, "row {i}: {} vs {}", g.1, w.1);
            }
        }
        // The coarse "constant" is sqrt(aggregate size): P_tent B_c = 1.
        let b_c = &tent.coarse_null_space.values;
        assert!((b_c[0] - 3.0_f64.sqrt()).abs() < 1e-14);
        assert!((b_c[1] - 2.0_f64.sqrt()).abs() < 1e-14);
        assert!((b_c[2] - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_sa_elasticity_prolongators_reproduce_rigid_body_modes_on_every_level() {
        // Free-floating cube: A B = 0 exactly, so every smoothed prolongator
        // keeps the (coarse) rigid-body modes in its range.  Seeing the six
        // fine-level modes in the range of P_0 P_1 requires that the coarse
        // near-null space R of level 0 was carried to level 1.
        let (a, coords, _) = cube_elasticity(6, false);
        assert_eq!(a.nrows, 3 * coords.len());
        let sa = SmoothedAggregationAmg {
            coarse_cutoff: 80,
            ..SmoothedAggregationAmg::new()
        };
        let hier = sa.build(&a, Some(&coords));
        assert!(hier.levels.len() >= 3, "levels: {}", hier.levels.len());
        let modes = rigid_body_modes_3d(&coords, 3);
        assert_eq!(modes.len(), 6);
        let p0 = &hier.levels[0].p;
        let composite = spmm(p0, &hier.levels[1].p);
        for (m, mode) in modes.iter().enumerate() {
            let r0 = range_residual(p0, mode);
            assert!(r0 < 1e-9, "mode {m} not in range(P0): residual {r0}");
            let r01 = range_residual(&composite, mode);
            assert!(r01 < 1e-9, "mode {m} not in range(P0 P1): residual {r01}");
        }
    }

    #[test]
    fn test_sa_elasticity_path_is_taken_with_six_dofs_per_coarse_node() {
        let (a, coords, _) = cube_elasticity(5, false);
        let sa = SmoothedAggregationAmg::new();
        let hier = sa.build(&a, Some(&coords));
        assert!(hier.levels.len() >= 2);
        // Recompute the first-level aggregation independently.
        let node_offsets: Vec<usize> = (0..=coords.len()).map(|i| 3 * i).collect();
        let dof_to_node = dof_to_node_map(&node_offsets);
        let strong = node_strong_connections(&a, &node_offsets, &dof_to_node, sa.sa_theta);
        let aggregate = aggregate_connected_nodes(&strong, coords.len());
        let n_aggs = aggregate.iter().max().map_or(0, |m| m + 1);
        let n_coarse = hier.levels[1].a.nrows;
        // The scalar path would produce one coarse DOF per (DOF) aggregate.
        assert_eq!(n_coarse, hier.levels[0].p.ncols);
        assert!(n_aggs >= 4 && n_aggs * 3 <= n_coarse && n_coarse <= n_aggs * 6);
        assert!(
            n_coarse * 10 >= n_aggs * 50,
            "most aggregates should carry all six rigid-body modes: {n_coarse} coarse DOFs for {n_aggs} aggregates"
        );
        assert!(
            n_coarse < a.nrows / 2,
            "level 1 has {n_coarse} of {} DOFs",
            a.nrows
        );
    }

    #[test]
    fn test_sa_elasticity_dirichlet_nodes_get_empty_prolongator_rows() {
        let (a, coords, fixed) = cube_elasticity(5, true);
        let sa = SmoothedAggregationAmg::new();
        let hier = sa.build(&a, Some(&coords));
        assert!(hier.levels.len() >= 2);
        let p = &hier.levels[0].p;
        assert!(!fixed.is_empty());
        for &dof in &fixed {
            assert!(
                dense_row(p, dof).iter().all(|e| e.1 == 0.0),
                "clamped DOF {dof} must not be interpolated"
            );
        }
        // Every free DOF is interpolated, and the coarse operators are SPD-like.
        let free_rows_with_entries = (0..p.nrows)
            .filter(|i| !fixed.contains(i))
            .filter(|&i| p.row_offsets[i + 1] > p.row_offsets[i])
            .count();
        assert_eq!(free_rows_with_entries, p.nrows - fixed.len());
        for level in &hier.levels {
            assert!(level.diag.iter().all(|d| *d > 0.0 && d.is_finite()));
        }
    }

    #[test]
    fn test_sa_elasticity_pcg_converges_in_few_iterations() {
        let (a, coords, _) = cube_elasticity(6, true);
        let n = a.nrows;
        // Load: unit force in x on the top face nodes.
        let mut b = vec![0.0_f64; n];
        for (node, c) in coords.iter().enumerate() {
            if c[2] == 1.0 {
                b[3 * node] = 1.0;
                b[3 * node + 2] = -0.5;
            }
        }
        let sa = SmoothedAggregationAmg::new();
        let hierarchy = sa.build(&a, Some(&coords));
        assert!(hierarchy.levels.len() >= 2);
        let precond = AmgPreconditioner {
            hierarchy,
            cycle_kind: crate::solvers::amg::cycle::CycleKind::V,
            pcg: ParallelPcgSolver::new(500, 1e-10),
        };
        let solver = PcgWithPrecond::new(precond, 100, 1e-8);
        let mut x = vec![0.0_f64; n];
        let stats = solver.solve(&a, &b, &mut x);
        assert!(stats.converged, "PCG+SA-AMG did not converge: {stats:?}");
        assert!(stats.iterations < 50, "iterations = {}", stats.iterations);
        // True residual check.
        let mut ax = vec![0.0_f64; n];
        a.spmv(&x, &mut ax);
        let res: f64 = ax
            .iter()
            .zip(b.iter())
            .map(|(p, q)| (p - q).powi(2))
            .sum::<f64>()
            .sqrt();
        let b_norm: f64 = b.iter().map(|v| v * v).sum::<f64>().sqrt();
        assert!(res / b_norm < 1e-7, "relative residual {}", res / b_norm);
        // A plain Jacobi-PCG needs many more iterations.
        let jacobi = ParallelPcgSolver::new(2000, 1e-8);
        let mut x_j = vec![0.0_f64; n];
        let jac_stats = jacobi.solve(&a, &b, &mut x_j);
        assert!(
            stats.iterations * 3 < jac_stats.iterations,
            "SA-AMG: {}, Jacobi: {}",
            stats.iterations,
            jac_stats.iterations
        );
    }

    #[test]
    fn test_sa_mismatched_coordinates_fall_back_to_the_scalar_path() {
        let (a, coords, _) = cube_elasticity(4, true);
        let sa = SmoothedAggregationAmg::new();
        let scalar = sa.build(&a, None);
        // One coordinate too few: 3 * len != nrows, so no elasticity path (and
        // no out-of-bounds indexing, which the old code hit).
        let short = &coords[..coords.len() - 1];
        let fallback = sa.build(&a, Some(short));
        assert_eq!(scalar.levels.len(), fallback.levels.len());
        for (l1, l2) in scalar.levels.iter().zip(fallback.levels.iter()) {
            assert_eq!(l1.a.nrows, l2.a.nrows);
        }
        let empty: Vec<[f64; 3]> = Vec::new();
        assert_eq!(sa.build(&a, Some(&empty)).levels.len(), scalar.levels.len());
    }

    #[test]
    fn test_sa_jacobi_omega_zero_gives_the_orthonormal_tentative_prolongator() {
        let a = make_1d_poisson(40);
        let sa = SmoothedAggregationAmg {
            coarse_cutoff: 8,
            jacobi_omega: 0.0,
            ..SmoothedAggregationAmg::new()
        };
        let hier = sa.build(&a, None);
        let p = &hier.levels[0].p;
        let nc = p.ncols;
        for c1 in 0..nc {
            for c2 in 0..nc {
                let mut dot = 0.0;
                for i in 0..p.nrows {
                    let row = dense_row(p, i);
                    let v1 = row.iter().find(|e| e.0 == c1).map_or(0.0, |e| e.1);
                    let v2 = row.iter().find(|e| e.0 == c2).map_or(0.0, |e| e.1);
                    dot += v1 * v2;
                }
                let want = if c1 == c2 { 1.0 } else { 0.0 };
                assert!((dot - want).abs() < 1e-13, "PtP[{c1}][{c2}] = {dot}");
            }
        }
    }
}
