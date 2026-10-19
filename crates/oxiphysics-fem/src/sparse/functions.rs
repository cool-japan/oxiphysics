//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::collections::HashMap;

use super::types::{CsrMatrix, QuadTreeNode};

/// Compute an Approximate Minimum Degree (AMD) ordering permutation for a
/// symmetric sparse matrix.
///
/// Returns a permutation vector `perm` such that `perm[i]` is the original
/// index of the `i`-th node in the reordered system.  The reordering tends
/// to reduce fill-in during LU/Cholesky factorization.
///
/// This implements the greedy minimum-degree heuristic: at each step, select
/// the node with the fewest structural neighbours (minimum degree) and
/// eliminate it.
pub fn amd_ordering(a: &CsrMatrix) -> Vec<usize> {
    let n = a.nrows;
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for row in 0..n {
        let start = a.row_ptr[row];
        let end = a.row_ptr[row + 1];
        for idx in start..end {
            let col = a.col_indices[idx];
            if col != row {
                adj[row].push(col);
                adj[col].push(row);
            }
        }
    }
    for v in adj.iter_mut() {
        v.sort_unstable();
        v.dedup();
    }
    let mut eliminated = vec![false; n];
    let mut perm = Vec::with_capacity(n);
    for _ in 0..n {
        let next = (0..n)
            .filter(|&v| !eliminated[v])
            .min_by_key(|&v| adj[v].iter().filter(|&&nb| !eliminated[nb]).count())
            .expect("at least one uneliminated node");
        perm.push(next);
        eliminated[next] = true;
        let live_neighbours: Vec<usize> = adj[next]
            .iter()
            .filter(|&&nb| !eliminated[nb])
            .copied()
            .collect();
        for &u in &live_neighbours {
            for &v in &live_neighbours {
                if u != v && !adj[u].contains(&v) {
                    adj[u].push(v);
                }
            }
        }
    }
    perm
}
/// Symmetrically permute a CSR matrix using permutation `perm`.
///
/// Returns `B` where `B[i,j] = A[perm[i\], perm[j]]`.
pub fn permute_matrix(a: &CsrMatrix, perm: &[usize]) -> CsrMatrix {
    let n = a.nrows;
    assert_eq!(
        n, a.ncols,
        "matrix must be square for symmetric permutation"
    );
    assert_eq!(perm.len(), n);
    let mut inv_perm = vec![0usize; n];
    for (new_idx, &old_idx) in perm.iter().enumerate() {
        inv_perm[old_idx] = new_idx;
    }
    let mut triplets = Vec::with_capacity(a.nnz());
    for (row, &perm_row) in perm.iter().enumerate() {
        let old_row = perm_row;
        let start = a.row_ptr[old_row];
        let end = a.row_ptr[old_row + 1];
        for idx in start..end {
            let old_col = a.col_indices[idx];
            let new_col = inv_perm[old_col];
            triplets.push((row, new_col, a.values[idx]));
        }
    }
    CsrMatrix::from_triplets(n, n, &triplets)
}
/// Apply a permutation to a vector: `y[i] = x[perm[i\]]`.
pub fn apply_permutation(x: &[f64], perm: &[usize]) -> Vec<f64> {
    perm.iter().map(|&i| x[i]).collect()
}
/// Apply the inverse permutation: `y[perm[i\]] = x[i]`.
pub fn apply_inverse_permutation(x: &[f64], perm: &[usize]) -> Vec<f64> {
    let n = x.len();
    let mut y = vec![0.0; n];
    for (i, &p) in perm.iter().enumerate() {
        y[p] = x[i];
    }
    y
}
/// Coarsen a symmetric sparse matrix using a greedy maximal independent set
/// (C/F splitting) heuristic.
///
/// Returns the set of coarse-grid node indices (C-points).  The remaining
/// nodes are fine-grid (F-points) and will be interpolated.
///
/// # Strategy
/// The function uses a strength-of-connection threshold: node `j` is
/// "strongly connected" to node `i` if `|a_ij| >= theta * max_k |a_ik|`.
/// A greedy MIS on the strong-connection graph defines the C-points.
pub fn amg_coarsen(a: &CsrMatrix) -> Vec<usize> {
    let n = a.nrows;
    let theta = 0.25;
    let mut max_off: Vec<f64> = vec![0.0; n];
    for (row, max_off_row) in max_off.iter_mut().enumerate() {
        let start = a.row_ptr[row];
        let end = a.row_ptr[row + 1];
        for idx in start..end {
            if a.col_indices[idx] != row {
                let v = a.values[idx].abs();
                if v > *max_off_row {
                    *max_off_row = v;
                }
            }
        }
    }
    let mut strong: Vec<Vec<usize>> = vec![Vec::new(); n];
    for row in 0..n {
        let start = a.row_ptr[row];
        let end = a.row_ptr[row + 1];
        let threshold = theta * max_off[row];
        for idx in start..end {
            let col = a.col_indices[idx];
            if col != row && a.values[idx].abs() >= threshold {
                strong[row].push(col);
            }
        }
    }
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| strong[b].len().cmp(&strong[a].len()));
    let mut coarse = Vec::new();
    let mut assigned = vec![false; n];
    for &node in &order {
        if assigned[node] {
            continue;
        }
        coarse.push(node);
        assigned[node] = true;
        for &nb in &strong[node] {
            assigned[nb] = true;
        }
    }
    for (i, &is_assigned) in assigned.iter().enumerate() {
        if !is_assigned {
            coarse.push(i);
        }
    }
    coarse.sort_unstable();
    coarse
}
/// Build a piecewise-constant prolongation operator P (fine × coarse).
///
/// Each fine-grid F-point is interpolated from the nearest C-point based on
/// the strength of connection.  C-points map to themselves with weight 1.
pub fn amg_prolongation(a: &CsrMatrix, coarse: &[usize]) -> CsrMatrix {
    let n_fine = a.nrows;
    let n_coarse = coarse.len();
    let mut coarse_map: Vec<Option<usize>> = vec![None; n_fine];
    for (ci, &fi) in coarse.iter().enumerate() {
        coarse_map[fi] = Some(ci);
    }
    let theta = 0.25;
    let mut max_off: Vec<f64> = vec![0.0; n_fine];
    for (row, max_off_row) in max_off.iter_mut().enumerate() {
        let start = a.row_ptr[row];
        let end = a.row_ptr[row + 1];
        for idx in start..end {
            if a.col_indices[idx] != row {
                let v = a.values[idx].abs();
                if v > *max_off_row {
                    *max_off_row = v;
                }
            }
        }
    }
    let mut triplets = Vec::new();
    for fine_row in 0..n_fine {
        if let Some(ci) = coarse_map[fine_row] {
            triplets.push((fine_row, ci, 1.0));
        } else {
            let start = a.row_ptr[fine_row];
            let end = a.row_ptr[fine_row + 1];
            let threshold = theta * max_off[fine_row];
            let mut den = 0.0f64;
            let mut c_contributions: Vec<(usize, f64)> = Vec::new();
            for idx in start..end {
                let col = a.col_indices[idx];
                if col != fine_row
                    && a.values[idx].abs() >= threshold
                    && let Some(coarse_idx) = coarse_map[col]
                {
                    let w = a.values[idx].abs();
                    c_contributions.push((coarse_idx, w));
                    den += w;
                }
            }
            if den > 1e-30 {
                for (ci, w) in c_contributions {
                    triplets.push((fine_row, ci, w / den));
                }
            } else {
                triplets.push((fine_row, 0, 1.0));
            }
        }
    }
    CsrMatrix::from_triplets(n_fine, n_coarse, &triplets)
}
/// Compute the Galerkin coarse-grid operator: A_c = P^T A P.
///
/// This is the standard Galerkin coarsening used in AMG.
pub fn amg_galerkin_coarse(a: &CsrMatrix, p: &CsrMatrix) -> CsrMatrix {
    let r = p.transpose();
    let ap = csr_mat_mul(a, p);
    csr_mat_mul(&r, &ap)
}
/// Sparse matrix-matrix multiplication: C = A * B.
///
/// Both matrices must be compatible: A is (m×k) and B is (k×n).
pub fn csr_mat_mul(a: &CsrMatrix, b: &CsrMatrix) -> CsrMatrix {
    assert_eq!(a.ncols, b.nrows, "inner dimensions must match");
    let m = a.nrows;
    let n = b.ncols;
    let mut triplets: Vec<(usize, usize, f64)> = Vec::new();
    for i in 0..m {
        let mut row_map: HashMap<usize, f64> = HashMap::new();
        let a_start = a.row_ptr[i];
        let a_end = a.row_ptr[i + 1];
        for a_idx in a_start..a_end {
            let k = a.col_indices[a_idx];
            let a_val = a.values[a_idx];
            let b_start = b.row_ptr[k];
            let b_end = b.row_ptr[k + 1];
            for b_idx in b_start..b_end {
                let j = b.col_indices[b_idx];
                *row_map.entry(j).or_insert(0.0) += a_val * b.values[b_idx];
            }
        }
        for (j, v) in row_map {
            if v.abs() > 1e-30 {
                triplets.push((i, j, v));
            }
        }
    }
    CsrMatrix::from_triplets(m, n, &triplets)
}
/// Perform one AMG two-level V-cycle to solve `A x = b`.
///
/// This is a one-shot helper: the coarsening, interpolation and Galerkin
/// operators are rebuilt on every call and the coarse system is only relaxed
/// by forward Gauss-Seidel sweeps, so the map `b ↦ x` is **not** a symmetric
/// operator.  To precondition a conjugate-gradient iteration use
/// [`crate::solvers::SymmetricAmg`], which builds the hierarchy once and applies
/// a symmetric positive definite V-cycle.
///
/// # Arguments
/// * `a`    – fine-grid operator
/// * `b`    – right-hand side
/// * `x0`  – initial guess
/// * `nu1`  – number of pre-smoothing steps (Jacobi)
/// * `nu2`  – number of post-smoothing steps (Jacobi)
///
/// Returns the updated solution vector.
pub fn amg_v_cycle(a: &CsrMatrix, b: &[f64], x0: &[f64], nu1: usize, nu2: usize) -> Vec<f64> {
    let n = a.nrows;
    let coarse = amg_coarsen(a);
    let p = amg_prolongation(a, &coarse);
    let ac = amg_galerkin_coarse(a, &p);
    let omega_j = 2.0 / 3.0;
    let mut x = x0.to_vec();
    for _ in 0..nu1 {
        x = jacobi_smooth(a, b, &x, omega_j);
    }
    let ax = a.mul_vec(&x);
    let r: Vec<f64> = (0..n).map(|i| b[i] - ax[i]).collect();
    let r_c = restrict_vector(&p, &r);
    let n_c = coarse.len();
    let e_c = gauss_seidel_solve(&ac, &r_c, &vec![0.0; n_c], 20);
    let e = prolongate_vector(&p, &e_c);
    for i in 0..n {
        x[i] += e[i];
    }
    for _ in 0..nu2 {
        x = jacobi_smooth(a, b, &x, omega_j);
    }
    x
}
/// Damped Jacobi smoother: one step.
pub(super) fn jacobi_smooth(a: &CsrMatrix, b: &[f64], x: &[f64], omega: f64) -> Vec<f64> {
    let n = a.nrows;
    let mut x_new = x.to_vec();
    for i in 0..n {
        let diag = a.diagonal(i);
        if diag.abs() < 1e-30 {
            continue;
        }
        let ax_i: f64 = {
            let start = a.row_ptr[i];
            let end = a.row_ptr[i + 1];
            let mut s = 0.0;
            for idx in start..end {
                if a.col_indices[idx] != i {
                    s += a.values[idx] * x[a.col_indices[idx]];
                }
            }
            s
        };
        x_new[i] = (1.0 - omega) * x[i] + omega * (b[i] - ax_i) / diag;
    }
    x_new
}
/// Restrict a fine-grid vector to the coarse grid: r_c = P^T r.
pub(super) fn restrict_vector(p: &CsrMatrix, r: &[f64]) -> Vec<f64> {
    let n_coarse = p.ncols;
    let mut r_c = vec![0.0f64; n_coarse];
    for (i, &r_i) in r.iter().enumerate().take(p.nrows) {
        let start = p.row_ptr[i];
        let end = p.row_ptr[i + 1];
        for idx in start..end {
            let j = p.col_indices[idx];
            r_c[j] += p.values[idx] * r_i;
        }
    }
    r_c
}
/// Prolongate a coarse-grid vector to the fine grid: e = P * e_c.
pub(super) fn prolongate_vector(p: &CsrMatrix, e_c: &[f64]) -> Vec<f64> {
    p.mul_vec(e_c)
}
/// Simple Gauss-Seidel solver for the coarse-grid system.
pub(super) fn gauss_seidel_solve(
    a: &CsrMatrix,
    b: &[f64],
    x0: &[f64],
    max_iter: usize,
) -> Vec<f64> {
    let n = a.nrows;
    let mut x = x0.to_vec();
    for _ in 0..max_iter {
        for i in 0..n {
            let diag = a.diagonal(i);
            if diag.abs() < 1e-30 {
                continue;
            }
            let start = a.row_ptr[i];
            let end = a.row_ptr[i + 1];
            let mut s = 0.0;
            for idx in start..end {
                let j = a.col_indices[idx];
                if j != i {
                    s += a.values[idx] * x[j];
                }
            }
            x[i] = (b[i] - s) / diag;
        }
    }
    x
}
/// Adaptively refine a quadtree based on element error indicators.
///
/// Nodes whose `error_indicator` exceeds `threshold` are refined, up to a
/// maximum refinement level `max_level`.
pub fn adaptive_refine_quadtree(node: &mut QuadTreeNode, threshold: f64, max_level: u32) {
    if node.level >= max_level {
        return;
    }
    if node.is_leaf() {
        if node.error_indicator >= threshold {
            node.refine();
        }
    } else if let Some(children) = node.children.as_mut() {
        for child in children.iter_mut() {
            adaptive_refine_quadtree(child, threshold, max_level);
        }
    }
}
/// Coarsen a quadtree: if all four siblings are leaves with error below
/// `threshold`, replace the parent with a single leaf.
///
/// Returns `true` if this node was coarsened.
pub fn coarsen_quadtree(node: &mut QuadTreeNode, threshold: f64) -> bool {
    if node.is_leaf() {
        return false;
    }
    if let Some(children) = node.children.as_mut() {
        for child in children.iter_mut() {
            coarsen_quadtree(child, threshold);
        }
    }
    let all_leaf_low = node.children.as_ref().is_some_and(|ch| {
        ch.iter()
            .all(|c| c.is_leaf() && c.error_indicator < threshold)
    });
    if all_leaf_low {
        node.children = None;
        return true;
    }
    false
}
/// Compute the maximum refinement level in a quadtree.
pub fn max_refinement_level(node: &QuadTreeNode) -> u32 {
    match &node.children {
        None => node.level,
        Some(ch) => ch
            .iter()
            .map(max_refinement_level)
            .max()
            .unwrap_or(node.level),
    }
}
/// Compute the 1-norm of a CSR matrix: max column sum of absolute values.
///
/// `||A||_1 = max_j sum_i |a_ij|`
pub fn csr_norm_1(a: &CsrMatrix) -> f64 {
    let mut col_sums = vec![0.0f64; a.ncols];
    for row in 0..a.nrows {
        let start = a.row_ptr[row];
        let end = a.row_ptr[row + 1];
        for idx in start..end {
            col_sums[a.col_indices[idx]] += a.values[idx].abs();
        }
    }
    col_sums.into_iter().fold(0.0_f64, f64::max)
}
/// Compute the infinity-norm of a CSR matrix: max row sum of absolute values.
///
/// `||A||_inf = max_i sum_j |a_ij|`
pub fn csr_norm_inf(a: &CsrMatrix) -> f64 {
    (0..a.nrows)
        .map(|row| {
            let start = a.row_ptr[row];
            let end = a.row_ptr[row + 1];
            a.values[start..end].iter().map(|v| v.abs()).sum::<f64>()
        })
        .fold(0.0_f64, f64::max)
}
/// Reverse Cuthill-McKee (RCM) ordering for bandwidth reduction.
///
/// Returns a permutation `perm` such that the reordered matrix `A[perm, perm]`
/// has a reduced bandwidth.  The algorithm performs a BFS from a peripheral
/// node, then reverses the resulting ordering.
///
/// # References
/// E. Cuthill and J. McKee, "Reducing the Bandwidth of Sparse Symmetric
/// Matrices", ACM 1969.
pub fn reverse_cuthill_mckee(a: &CsrMatrix) -> Vec<usize> {
    assert_eq!(a.nrows, a.ncols, "RCM requires a square symmetric matrix");
    let n = a.nrows;
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for row in 0..n {
        let start = a.row_ptr[row];
        let end = a.row_ptr[row + 1];
        for idx in start..end {
            let col = a.col_indices[idx];
            if col != row {
                adj[row].push(col);
                adj[col].push(row);
            }
        }
    }
    for v in adj.iter_mut() {
        v.sort_unstable();
        v.dedup();
    }
    let start_node = (0..n).min_by_key(|&i| adj[i].len()).unwrap_or(0);
    let mut visited = vec![false; n];
    let mut queue = std::collections::VecDeque::new();
    let mut order = Vec::with_capacity(n);
    queue.push_back(start_node);
    visited[start_node] = true;
    while let Some(node) = queue.pop_front() {
        order.push(node);
        let mut neighbours: Vec<usize> = adj[node]
            .iter()
            .filter(|&&nb| !visited[nb])
            .copied()
            .collect();
        neighbours.sort_by_key(|&nb| adj[nb].len());
        for nb in neighbours {
            if !visited[nb] {
                visited[nb] = true;
                queue.push_back(nb);
            }
        }
    }
    for (i, &v) in visited.iter().enumerate() {
        if !v {
            order.push(i);
        }
    }
    order.reverse();
    order
}
/// Compute the bandwidth of a sparse matrix under a given permutation.
///
/// Returns `max_{(i,j) in A} |new_i - new_j|` where `new_k = inv_perm[k]`.
pub fn bandwidth(a: &CsrMatrix, perm: &[usize]) -> usize {
    let n = a.nrows;
    let mut inv_perm = vec![0usize; n];
    for (new, &old) in perm.iter().enumerate() {
        inv_perm[old] = new;
    }
    let mut bw = 0usize;
    for row in 0..n {
        let start = a.row_ptr[row];
        let end = a.row_ptr[row + 1];
        for idx in start..end {
            let col = a.col_indices[idx];
            let new_row = inv_perm[row];
            let new_col = inv_perm[col];
            let diff = new_row.abs_diff(new_col);
            if diff > bw {
                bw = diff;
            }
        }
    }
    bw
}
/// Convert a full symmetric CSR matrix to upper-triangle-only storage.
///
/// Returns a new CSR matrix containing only the entries `(i, j)` with `j >= i`.
/// The diagonal and strictly upper triangle are retained.
pub fn to_upper_triangular(a: &CsrMatrix) -> CsrMatrix {
    let mut triplets = Vec::new();
    for row in 0..a.nrows {
        let start = a.row_ptr[row];
        let end = a.row_ptr[row + 1];
        for idx in start..end {
            let col = a.col_indices[idx];
            if col >= row {
                triplets.push((row, col, a.values[idx]));
            }
        }
    }
    CsrMatrix::from_triplets(a.nrows, a.ncols, &triplets)
}
/// Convert an upper-triangle-only CSR matrix back to full symmetric storage.
///
/// Each off-diagonal entry `(i, j)` in the upper triangle generates a
/// symmetric entry `(j, i)` with the same value.
pub fn from_upper_triangular(a: &CsrMatrix) -> CsrMatrix {
    let mut triplets = Vec::new();
    for row in 0..a.nrows {
        let start = a.row_ptr[row];
        let end = a.row_ptr[row + 1];
        for idx in start..end {
            let col = a.col_indices[idx];
            let val = a.values[idx];
            triplets.push((row, col, val));
            if col != row {
                triplets.push((col, row, val));
            }
        }
    }
    CsrMatrix::from_triplets(a.nrows, a.ncols, &triplets)
}
/// Check if a CSR matrix is numerically symmetric within a given tolerance.
///
/// Returns `true` if `|A[i,j] - A[j,i]| <= tol * max(|A[i,j]|, |A[j,i]|, 1)`
/// for all stored entries.
pub fn is_symmetric(a: &CsrMatrix, tol: f64) -> bool {
    for row in 0..a.nrows {
        let start = a.row_ptr[row];
        let end = a.row_ptr[row + 1];
        for idx in start..end {
            let col = a.col_indices[idx];
            let aij = a.values[idx];
            let aji = a.get(col, row);
            let scale = aij.abs().max(aji.abs()).max(1.0);
            if (aij - aji).abs() / scale > tol {
                return false;
            }
        }
    }
    true
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::sparse::*;
    #[test]
    fn test_csr_from_triplets() {
        let triplets = vec![
            (0, 0, 4.0),
            (0, 1, 1.0),
            (1, 0, 1.0),
            (1, 1, 3.0),
            (1, 1, 1.0),
        ];
        let m = CsrMatrix::from_triplets(2, 2, &triplets);
        assert_eq!(m.get(0, 0), 4.0);
        assert_eq!(m.get(0, 1), 1.0);
        assert_eq!(m.get(1, 0), 1.0);
        assert_eq!(m.get(1, 1), 4.0);
        assert_eq!(m.nnz(), 4);
    }
    #[test]
    fn test_csr_mul_vec() {
        let triplets = vec![(0, 0, 2.0), (0, 1, 1.0), (1, 0, 1.0), (1, 1, 3.0)];
        let m = CsrMatrix::from_triplets(2, 2, &triplets);
        let result = m.mul_vec(&[1.0, 2.0]);
        assert!((result[0] - 4.0).abs() < 1e-12);
        assert!((result[1] - 7.0).abs() < 1e-12);
    }
    #[test]
    fn test_csr_add_to() {
        let mut m = CsrMatrix::new(3, 3);
        m.add_to(0, 0, 1.0);
        m.add_to(0, 0, 2.0);
        m.add_to(1, 2, 5.0);
        assert_eq!(m.get(0, 0), 3.0);
        assert_eq!(m.get(1, 2), 5.0);
        assert_eq!(m.get(2, 2), 0.0);
    }
    #[test]
    fn test_csr_set() {
        let mut m = CsrMatrix::new(2, 2);
        m.set(0, 0, 10.0);
        m.set(0, 1, 20.0);
        m.set(1, 1, 30.0);
        assert_eq!(m.get(0, 0), 10.0);
        assert_eq!(m.get(0, 1), 20.0);
        assert_eq!(m.get(1, 0), 0.0);
        assert_eq!(m.get(1, 1), 30.0);
        m.set(0, 0, 99.0);
        assert_eq!(m.get(0, 0), 99.0);
    }
    #[test]
    fn test_sparse_vector() {
        let a = SparseVector::from_vec(vec![1.0, 2.0, 3.0]);
        let b = SparseVector::from_vec(vec![4.0, 5.0, 6.0]);
        assert!((a.dot(&b) - 32.0).abs() < 1e-12);
        assert!((a.norm() - 14.0_f64.sqrt()).abs() < 1e-12);
    }
    #[test]
    fn test_csr_transpose() {
        let triplets = vec![(0, 0, 1.0), (0, 1, 2.0), (1, 0, 3.0), (1, 1, 4.0)];
        let m = CsrMatrix::from_triplets(2, 2, &triplets);
        let mt = m.transpose();
        assert_eq!(mt.get(0, 0), 1.0);
        assert_eq!(mt.get(0, 1), 3.0);
        assert_eq!(mt.get(1, 0), 2.0);
        assert_eq!(mt.get(1, 1), 4.0);
    }
    #[test]
    fn test_csr_transpose_rectangular() {
        let triplets = vec![(0, 0, 1.0), (0, 1, 2.0), (0, 2, 3.0)];
        let m = CsrMatrix::from_triplets(1, 3, &triplets);
        let mt = m.transpose();
        assert_eq!(mt.nrows, 3);
        assert_eq!(mt.ncols, 1);
        assert_eq!(mt.get(0, 0), 1.0);
        assert_eq!(mt.get(1, 0), 2.0);
        assert_eq!(mt.get(2, 0), 3.0);
    }
    #[test]
    fn test_csr_add() {
        let a = CsrMatrix::from_triplets(2, 2, &[(0, 0, 1.0), (1, 1, 2.0)]);
        let b = CsrMatrix::from_triplets(2, 2, &[(0, 0, 10.0), (0, 1, 5.0)]);
        let c = a.add(&b);
        assert_eq!(c.get(0, 0), 11.0);
        assert_eq!(c.get(0, 1), 5.0);
        assert_eq!(c.get(1, 0), 0.0);
        assert_eq!(c.get(1, 1), 2.0);
    }
    #[test]
    fn test_csr_scale() {
        let mut m = CsrMatrix::from_triplets(2, 2, &[(0, 0, 3.0), (1, 1, 4.0)]);
        m.scale(2.0);
        assert_eq!(m.get(0, 0), 6.0);
        assert_eq!(m.get(1, 1), 8.0);
    }
    #[test]
    fn test_csr_scaled() {
        let m = CsrMatrix::from_triplets(2, 2, &[(0, 0, 3.0), (1, 1, 4.0)]);
        let ms = m.scaled(0.5);
        assert_eq!(ms.get(0, 0), 1.5);
        assert_eq!(ms.get(1, 1), 2.0);
        assert_eq!(m.get(0, 0), 3.0);
    }
    #[test]
    fn test_mul_vec_axpby() {
        let m =
            CsrMatrix::from_triplets(2, 2, &[(0, 0, 2.0), (0, 1, 1.0), (1, 0, 1.0), (1, 1, 3.0)]);
        let x = vec![1.0, 2.0];
        let mut y = vec![10.0, 20.0];
        m.mul_vec_axpby(&x, &mut y, 2.0, 0.5);
        assert!((y[0] - 13.0).abs() < 1e-12);
        assert!((y[1] - 24.0).abs() < 1e-12);
    }
    #[test]
    fn test_symmetric_mul_vec() {
        let m = CsrMatrix::from_triplets(2, 2, &[(0, 0, 4.0), (0, 1, 2.0), (1, 1, 3.0)]);
        let x = vec![1.0, 2.0];
        let y = m.symmetric_mul_vec(&x);
        assert!((y[0] - 8.0).abs() < 1e-12, "y[0] = {}", y[0]);
        assert!((y[1] - 8.0).abs() < 1e-12, "y[1] = {}", y[1]);
    }
    #[test]
    fn test_diagonal_vec() {
        let m = CsrMatrix::from_triplets(
            3,
            3,
            &[
                (0, 0, 5.0),
                (0, 1, 1.0),
                (1, 1, 7.0),
                (2, 0, 2.0),
                (2, 2, 9.0),
            ],
        );
        let diag = m.diagonal_vec();
        assert_eq!(diag, vec![5.0, 7.0, 9.0]);
    }
    #[test]
    fn test_frobenius_norm() {
        let m = CsrMatrix::from_triplets(1, 2, &[(0, 0, 3.0), (0, 1, 4.0)]);
        assert!((m.frobenius_norm() - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_csc_from_triplets() {
        let triplets = vec![(0, 0, 1.0), (0, 1, 2.0), (1, 0, 3.0), (1, 1, 4.0)];
        let m = CscMatrix::from_triplets(2, 2, &triplets);
        assert_eq!(m.get(0, 0), 1.0);
        assert_eq!(m.get(0, 1), 2.0);
        assert_eq!(m.get(1, 0), 3.0);
        assert_eq!(m.get(1, 1), 4.0);
        assert_eq!(m.nnz(), 4);
    }
    #[test]
    fn test_csc_mul_vec() {
        let triplets = vec![(0, 0, 2.0), (0, 1, 1.0), (1, 0, 1.0), (1, 1, 3.0)];
        let m = CscMatrix::from_triplets(2, 2, &triplets);
        let y = m.mul_vec(&[1.0, 2.0]);
        assert!((y[0] - 4.0).abs() < 1e-12);
        assert!((y[1] - 7.0).abs() < 1e-12);
    }
    #[test]
    fn test_csc_transpose() {
        let triplets = vec![(0, 0, 1.0), (0, 1, 2.0), (1, 0, 3.0)];
        let m = CscMatrix::from_triplets(2, 2, &triplets);
        let mt = m.transpose();
        assert_eq!(mt.get(0, 0), 1.0);
        assert_eq!(mt.get(1, 0), 2.0);
        assert_eq!(mt.get(0, 1), 3.0);
    }
    #[test]
    fn test_csr_to_csc_round_trip() {
        let triplets = vec![(0, 0, 1.0), (0, 1, 2.0), (1, 0, 3.0), (1, 1, 4.0)];
        let csr = CsrMatrix::from_triplets(2, 2, &triplets);
        let csc = csr.to_csc();
        let csr2 = csc.to_csr();
        for r in 0..2 {
            for c in 0..2 {
                assert!(
                    (csr.get(r, c) - csr2.get(r, c)).abs() < 1e-12,
                    "mismatch at ({r},{c}): {} vs {}",
                    csr.get(r, c),
                    csr2.get(r, c)
                );
            }
        }
    }
    #[test]
    fn test_csc_duplicate_sum() {
        let triplets = vec![(0, 0, 1.0), (0, 0, 2.0)];
        let m = CscMatrix::from_triplets(1, 1, &triplets);
        assert_eq!(m.get(0, 0), 3.0);
    }
    #[test]
    fn test_ilu0_diagonal_matrix() {
        let m = CsrMatrix::from_triplets(3, 3, &[(0, 0, 4.0), (1, 1, 5.0), (2, 2, 6.0)]);
        let ilu = Ilu0Preconditioner::new(&m);
        let rhs = vec![8.0, 15.0, 18.0];
        let z = ilu.solve(&rhs);
        assert!((z[0] - 2.0).abs() < 1e-10, "z[0] = {}", z[0]);
        assert!((z[1] - 3.0).abs() < 1e-10, "z[1] = {}", z[1]);
        assert!((z[2] - 3.0).abs() < 1e-10, "z[2] = {}", z[2]);
    }
    #[test]
    fn test_ilu0_tridiagonal() {
        let triplets = vec![
            (0, 0, 4.0),
            (0, 1, -1.0),
            (1, 0, -1.0),
            (1, 1, 4.0),
            (1, 2, -1.0),
            (2, 1, -1.0),
            (2, 2, 4.0),
        ];
        let m = CsrMatrix::from_triplets(3, 3, &triplets);
        let ilu = Ilu0Preconditioner::new(&m);
        let x_exact = vec![1.0, 2.0, 3.0];
        let rhs = m.mul_vec(&x_exact);
        let z = ilu.solve(&rhs);
        for i in 0..3 {
            assert!(
                (z[i] - x_exact[i]).abs() < 1e-10,
                "z[{i}] = {}, expected {}",
                z[i],
                x_exact[i]
            );
        }
    }
    #[test]
    fn test_ilu0_reduces_residual() {
        let triplets = vec![
            (0, 0, 10.0),
            (0, 1, 1.0),
            (0, 2, 0.5),
            (1, 0, 1.0),
            (1, 1, 8.0),
            (1, 2, 1.0),
            (2, 0, 0.5),
            (2, 1, 1.0),
            (2, 2, 6.0),
        ];
        let m = CsrMatrix::from_triplets(3, 3, &triplets);
        let ilu = Ilu0Preconditioner::new(&m);
        let rhs = vec![1.0, 2.0, 3.0];
        let z = ilu.solve(&rhs);
        let az = m.mul_vec(&z);
        let residual: f64 = rhs
            .iter()
            .zip(az.iter())
            .map(|(b, a)| (b - a) * (b - a))
            .sum::<f64>()
            .sqrt();
        let rhs_norm: f64 = rhs.iter().map(|b| b * b).sum::<f64>().sqrt();
        assert!(
            residual / rhs_norm < 0.1,
            "relative residual = {}",
            residual / rhs_norm
        );
    }
    #[test]
    fn test_sparse_vector_axpby() {
        let mut a = SparseVector::from_vec(vec![1.0, 2.0, 3.0]);
        let b = SparseVector::from_vec(vec![4.0, 5.0, 6.0]);
        a.axpby(2.0, &b, 3.0);
        assert!((a.data[0] - 14.0).abs() < 1e-12);
        assert!((a.data[1] - 19.0).abs() < 1e-12);
        assert!((a.data[2] - 24.0).abs() < 1e-12);
    }
    #[test]
    fn test_sparse_vector_scale() {
        let mut a = SparseVector::from_vec(vec![2.0, 4.0, 6.0]);
        a.scale(0.5);
        assert_eq!(a.data, vec![1.0, 2.0, 3.0]);
    }
    #[test]
    fn test_csr_add_identity() {
        let a = CsrMatrix::from_triplets(2, 2, &[(0, 0, 3.0), (1, 1, 7.0)]);
        let zero = CsrMatrix::from_triplets(2, 2, &[]);
        let c = a.add(&zero);
        assert_eq!(c.get(0, 0), 3.0);
        assert_eq!(c.get(1, 1), 7.0);
        assert_eq!(c.nnz(), 2);
    }
    #[test]
    fn test_csc_csr_mul_vec_consistency() {
        let triplets = vec![
            (0, 0, 2.0),
            (0, 1, 3.0),
            (1, 0, 4.0),
            (1, 1, 5.0),
            (2, 0, 6.0),
            (2, 1, 7.0),
        ];
        let csr = CsrMatrix::from_triplets(3, 2, &triplets);
        let csc = CscMatrix::from_triplets(3, 2, &triplets);
        let x = vec![1.0, 2.0];
        let y_csr = csr.mul_vec(&x);
        let y_csc = csc.mul_vec(&x);
        for i in 0..3 {
            assert!(
                (y_csr[i] - y_csc[i]).abs() < 1e-12,
                "mismatch at {i}: CSR={} CSC={}",
                y_csr[i],
                y_csc[i]
            );
        }
    }
    #[test]
    fn test_amd_identity_on_diagonal() {
        let m =
            CsrMatrix::from_triplets(4, 4, &[(0, 0, 1.0), (1, 1, 2.0), (2, 2, 3.0), (3, 3, 4.0)]);
        let perm = amd_ordering(&m);
        assert_eq!(perm.len(), 4);
        let mut sorted = perm.clone();
        sorted.sort_unstable();
        assert_eq!(sorted, vec![0, 1, 2, 3]);
    }
    #[test]
    fn test_amd_tridiagonal() {
        let triplets = vec![
            (0, 0, 4.0),
            (0, 1, -1.0),
            (1, 0, -1.0),
            (1, 1, 4.0),
            (1, 2, -1.0),
            (2, 1, -1.0),
            (2, 2, 4.0),
            (2, 3, -1.0),
            (3, 2, -1.0),
            (3, 3, 4.0),
        ];
        let m = CsrMatrix::from_triplets(4, 4, &triplets);
        let perm = amd_ordering(&m);
        assert_eq!(perm.len(), 4);
        let mut check = [false; 4];
        for &p in &perm {
            check[p] = true;
        }
        assert!(check.iter().all(|&v| v));
    }
    #[test]
    fn test_amd_permutation_matrix_apply() {
        let m = CsrMatrix::from_triplets(3, 3, &[(0, 0, 10.0), (1, 1, 20.0), (2, 2, 30.0)]);
        let perm = vec![2, 0, 1];
        let pm = permute_matrix(&m, &perm);
        assert!((pm.get(0, 0) - 30.0).abs() < 1e-12);
        assert!((pm.get(1, 1) - 10.0).abs() < 1e-12);
        assert!((pm.get(2, 2) - 20.0).abs() < 1e-12);
    }
    #[test]
    fn test_amd_permute_vector() {
        let perm = vec![2, 0, 1];
        let x = vec![10.0, 20.0, 30.0];
        let px = apply_permutation(&x, &perm);
        assert_eq!(px, vec![30.0, 10.0, 20.0]);
    }
    #[test]
    fn test_iluk_level0_matches_ilu0() {
        let triplets = vec![
            (0, 0, 10.0),
            (0, 1, 1.0),
            (0, 2, 0.5),
            (1, 0, 1.0),
            (1, 1, 8.0),
            (1, 2, 1.0),
            (2, 0, 0.5),
            (2, 1, 1.0),
            (2, 2, 6.0),
        ];
        let m = CsrMatrix::from_triplets(3, 3, &triplets);
        let iluk = IlukPreconditioner::new(&m, 0);
        let ilu0 = Ilu0Preconditioner::new(&m);
        let rhs = vec![1.0, 2.0, 3.0];
        let z_k = iluk.solve(&rhs);
        let z_0 = ilu0.solve(&rhs);
        for i in 0..3 {
            assert!(
                (z_k[i] - z_0[i]).abs() < 1e-10,
                "ILU(0) and ILU(k=0) differ at {i}: {} vs {}",
                z_k[i],
                z_0[i]
            );
        }
    }
    #[test]
    fn test_iluk_level1_diagonal() {
        let m = CsrMatrix::from_triplets(3, 3, &[(0, 0, 4.0), (1, 1, 5.0), (2, 2, 6.0)]);
        let iluk = IlukPreconditioner::new(&m, 1);
        let rhs = vec![8.0, 15.0, 18.0];
        let z = iluk.solve(&rhs);
        assert!((z[0] - 2.0).abs() < 1e-10);
        assert!((z[1] - 3.0).abs() < 1e-10);
        assert!((z[2] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_iluk_level1_reduces_residual() {
        let triplets = vec![
            (0, 0, 10.0),
            (0, 2, 1.0),
            (1, 1, 8.0),
            (1, 3, 1.0),
            (2, 0, 1.0),
            (2, 2, 6.0),
            (3, 1, 1.0),
            (3, 3, 5.0),
        ];
        let m = CsrMatrix::from_triplets(4, 4, &triplets);
        let iluk = IlukPreconditioner::new(&m, 1);
        let rhs = vec![1.0, 2.0, 1.0, 2.0];
        let z = iluk.solve(&rhs);
        let az = m.mul_vec(&z);
        let res: f64 = rhs
            .iter()
            .zip(az.iter())
            .map(|(b, a)| (b - a).powi(2))
            .sum::<f64>()
            .sqrt();
        let rhs_norm: f64 = rhs.iter().map(|b| b * b).sum::<f64>().sqrt();
        assert!(
            res / rhs_norm < 0.2,
            "ILU(1) relative residual = {}",
            res / rhs_norm
        );
    }
    #[test]
    fn test_amg_coarsen_returns_valid_coarse_set() {
        let triplets = vec![
            (0, 0, 4.0),
            (0, 1, -1.0),
            (1, 0, -1.0),
            (1, 1, 4.0),
            (1, 2, -1.0),
            (2, 1, -1.0),
            (2, 2, 4.0),
            (2, 3, -1.0),
            (3, 2, -1.0),
            (3, 3, 4.0),
        ];
        let m = CsrMatrix::from_triplets(4, 4, &triplets);
        let coarse = amg_coarsen(&m);
        assert!(!coarse.is_empty(), "Coarse set must not be empty");
        assert!(
            coarse.len() <= m.nrows,
            "Coarse set cannot exceed total nodes"
        );
    }
    #[test]
    fn test_amg_prolongation_dimensions() {
        let m = CsrMatrix::from_triplets(
            4,
            4,
            &[
                (0, 0, 4.0),
                (0, 1, -1.0),
                (1, 0, -1.0),
                (1, 1, 4.0),
                (1, 2, -1.0),
                (2, 1, -1.0),
                (2, 2, 4.0),
                (2, 3, -1.0),
                (3, 2, -1.0),
                (3, 3, 4.0),
            ],
        );
        let coarse = amg_coarsen(&m);
        let p = amg_prolongation(&m, &coarse);
        assert_eq!(p.nrows, m.nrows, "Prolongation rows = fine nodes");
        assert_eq!(p.ncols, coarse.len(), "Prolongation cols = coarse nodes");
    }
    #[test]
    fn test_amg_galerkin_coarse_operator() {
        let m = CsrMatrix::from_triplets(
            4,
            4,
            &[
                (0, 0, 4.0),
                (0, 1, -1.0),
                (1, 0, -1.0),
                (1, 1, 4.0),
                (1, 2, -1.0),
                (2, 1, -1.0),
                (2, 2, 4.0),
                (2, 3, -1.0),
                (3, 2, -1.0),
                (3, 3, 4.0),
            ],
        );
        let coarse = amg_coarsen(&m);
        let p = amg_prolongation(&m, &coarse);
        let ac = amg_galerkin_coarse(&m, &p);
        assert_eq!(ac.nrows, coarse.len());
        assert_eq!(ac.ncols, coarse.len());
    }
    #[test]
    fn test_amg_two_level_v_cycle_reduces_residual() {
        let n = 8;
        let mut triplets = vec![];
        for i in 0..n {
            triplets.push((i, i, 2.0));
            if i + 1 < n {
                triplets.push((i, i + 1, -1.0));
            }
            if i > 0 {
                triplets.push((i, i - 1, -1.0));
            }
        }
        let a = CsrMatrix::from_triplets(n, n, &triplets);
        let rhs = vec![1.0; n];
        let x0 = vec![0.0; n];
        let x = amg_v_cycle(&a, &rhs, &x0, 2, 2);
        let ax = a.mul_vec(&x);
        let res: f64 = rhs
            .iter()
            .zip(ax.iter())
            .map(|(b, a)| (b - a).powi(2))
            .sum::<f64>()
            .sqrt();
        let rhs_norm: f64 = rhs.iter().map(|b| b * b).sum::<f64>().sqrt();
        assert!(
            res / rhs_norm < 0.5,
            "AMG V-cycle relative residual = {}",
            res / rhs_norm
        );
    }
    #[test]
    fn test_amr_quad_tree_root() {
        let root = QuadTreeNode::new_root(0.0, 0.0, 1.0, 1.0);
        assert_eq!(root.level, 0);
        assert!(root.children.is_none());
        assert!(root.is_leaf());
    }
    #[test]
    fn test_amr_quad_tree_refine() {
        let mut root = QuadTreeNode::new_root(0.0, 0.0, 1.0, 1.0);
        root.refine();
        assert!(!root.is_leaf());
        let children = root.children.as_ref().unwrap();
        assert_eq!(children.len(), 4);
        for child in children.iter() {
            assert_eq!(child.level, 1);
        }
    }
    #[test]
    fn test_amr_quad_tree_leaf_count() {
        let mut root = QuadTreeNode::new_root(0.0, 0.0, 1.0, 1.0);
        assert_eq!(root.leaf_count(), 1);
        root.refine();
        assert_eq!(root.leaf_count(), 4);
    }
    #[test]
    fn test_amr_quad_tree_refine_nested() {
        let mut root = QuadTreeNode::new_root(0.0, 0.0, 1.0, 1.0);
        root.refine();
        if let Some(children) = root.children.as_mut() {
            children[0].refine();
        }
        assert_eq!(root.leaf_count(), 7);
    }
    #[test]
    fn test_amr_mesh_cell_size() {
        let root = QuadTreeNode::new_root(0.0, 0.0, 4.0, 4.0);
        assert!((root.cell_width() - 4.0).abs() < 1e-12);
        assert!((root.cell_height() - 4.0).abs() < 1e-12);
    }
    #[test]
    fn test_amr_mesh_child_positions() {
        let root = QuadTreeNode::new_root(0.0, 0.0, 2.0, 2.0);
        let mut r = root.clone();
        r.refine();
        let children = r.children.as_ref().unwrap();
        assert!((children[0].x0 - 0.0).abs() < 1e-12);
        assert!((children[0].y0 - 0.0).abs() < 1e-12);
        assert!((children[0].cell_width() - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_amr_error_indicator() {
        let mut root = QuadTreeNode::new_root(0.0, 0.0, 1.0, 1.0);
        root.error_indicator = 0.5;
        root.refine();
        let children = root.children.as_ref().unwrap();
        for c in children.iter() {
            assert!((c.error_indicator - 0.0).abs() < 1e-12);
        }
    }
    #[test]
    fn test_amr_adaptive_refine() {
        let mut root = QuadTreeNode::new_root(0.0, 0.0, 1.0, 1.0);
        root.error_indicator = 2.0;
        adaptive_refine_quadtree(&mut root, 1.0, 1);
        assert!(!root.is_leaf(), "High-error node should be refined");
    }
    #[test]
    fn test_amr_no_refine_below_threshold() {
        let mut root = QuadTreeNode::new_root(0.0, 0.0, 1.0, 1.0);
        root.error_indicator = 0.1;
        adaptive_refine_quadtree(&mut root, 1.0, 3);
        assert!(root.is_leaf(), "Low-error node should not be refined");
    }
    #[test]
    fn test_icc_diagonal_matrix() {
        let m = CsrMatrix::from_triplets(3, 3, &[(0, 0, 4.0), (1, 1, 9.0), (2, 2, 16.0)]);
        let icc = IccPreconditioner::new(&m);
        let rhs = vec![8.0, 27.0, 48.0];
        let z = icc.solve(&rhs);
        assert!((z[0] - 2.0).abs() < 1e-9, "z[0] = {}", z[0]);
        assert!((z[1] - 3.0).abs() < 1e-9, "z[1] = {}", z[1]);
        assert!((z[2] - 3.0).abs() < 1e-9, "z[2] = {}", z[2]);
    }
    #[test]
    fn test_icc_tridiagonal_spd() {
        let triplets = vec![
            (0, 0, 4.0),
            (0, 1, -1.0),
            (1, 0, -1.0),
            (1, 1, 4.0),
            (1, 2, -1.0),
            (2, 1, -1.0),
            (2, 2, 4.0),
        ];
        let m = CsrMatrix::from_triplets(3, 3, &triplets);
        let icc = IccPreconditioner::new(&m);
        let x_exact = vec![1.0, 1.0, 1.0];
        let rhs = m.mul_vec(&x_exact);
        let z = icc.solve(&rhs);
        let residual: f64 = {
            let mz = m.mul_vec(&z);
            rhs.iter()
                .zip(mz.iter())
                .map(|(b, a)| (b - a).powi(2))
                .sum::<f64>()
                .sqrt()
        };
        let rhs_norm: f64 = rhs.iter().map(|b| b * b).sum::<f64>().sqrt();
        assert!(
            residual / rhs_norm < 0.2,
            "ICC relative residual = {}",
            residual / rhs_norm
        );
    }
    #[test]
    fn test_icc_nnz_le_lower_triangle_nnz() {
        let triplets = vec![
            (0, 0, 5.0),
            (0, 1, -1.0),
            (0, 2, -1.0),
            (1, 0, -1.0),
            (1, 1, 5.0),
            (1, 2, -1.0),
            (2, 0, -1.0),
            (2, 1, -1.0),
            (2, 2, 5.0),
        ];
        let m = CsrMatrix::from_triplets(3, 3, &triplets);
        let icc = IccPreconditioner::new(&m);
        assert_eq!(icc.nnz(), 6, "ICC nnz should equal lower-triangle nnz");
    }
    /// Symmetric tridiagonal CSR matrix `tridiag(off, diag, off)` built directly
    /// from its arrays (no hashing), so large instances are cheap to create.
    fn tridiagonal_csr(n: usize, diag: f64, off: f64) -> CsrMatrix {
        let mut row_ptr = Vec::with_capacity(n + 1);
        let mut col_indices = Vec::with_capacity(3 * n);
        let mut values = Vec::with_capacity(3 * n);
        row_ptr.push(0);
        for i in 0..n {
            if i > 0 {
                col_indices.push(i - 1);
                values.push(off);
            }
            col_indices.push(i);
            values.push(diag);
            if i + 1 < n {
                col_indices.push(i + 1);
                values.push(off);
            }
            row_ptr.push(col_indices.len());
        }
        CsrMatrix {
            row_ptr,
            col_indices,
            values,
            nrows: n,
            ncols: n,
        }
    }
    /// 5-point Laplacian on an `m x m` grid (Dirichlet), `n = m * m`.
    fn poisson_2d_csr(m: usize) -> CsrMatrix {
        let mut triplets = Vec::new();
        for gy in 0..m {
            for gx in 0..m {
                let i = gy * m + gx;
                triplets.push((i, i, 4.0));
                if gx > 0 {
                    triplets.push((i, i - 1, -1.0));
                }
                if gx + 1 < m {
                    triplets.push((i, i + 1, -1.0));
                }
                if gy > 0 {
                    triplets.push((i, i - m, -1.0));
                }
                if gy + 1 < m {
                    triplets.push((i, i + m, -1.0));
                }
            }
        }
        CsrMatrix::from_triplets(m * m, m * m, &triplets)
    }
    /// Entry `(i, j)` of `L L^T` for the lower-triangular factor of `icc`.
    fn icc_product_entry(icc: &IccPreconditioner, i: usize, j: usize) -> f64 {
        let row = |r: usize| -> Vec<(usize, f64)> {
            (icc.row_ptr[r]..icc.row_ptr[r + 1])
                .map(|p| (icc.col_indices[p], icc.l_values[p]))
                .collect()
        };
        let row_i = row(i);
        let row_j = row(j);
        let mut sum = 0.0;
        for &(ci, vi) in &row_i {
            for &(cj, vj) in &row_j {
                if ci == cj {
                    sum += vi * vj;
                }
            }
        }
        sum
    }
    #[test]
    fn test_icc_reproduces_matrix_on_its_sparsity_pattern() {
        // Defining property of IC(0): (L L^T)_ij = A_ij wherever A_ij is stored.
        for a in [poisson_2d_csr(5), tridiagonal_csr(9, 3.0, -1.0)] {
            let icc = IccPreconditioner::new(&a);
            assert_eq!(icc.shift(), 0.0);
            assert!(!icc.used_diagonal_fallback());
            for i in 0..a.nrows {
                for p in a.row_ptr[i]..a.row_ptr[i + 1] {
                    let j = a.col_indices[p];
                    let reproduced = icc_product_entry(&icc, i, j);
                    assert!(
                        (reproduced - a.values[p]).abs() < 1e-12,
                        "(L L^T)[{i},{j}] = {reproduced}, A = {}",
                        a.values[p]
                    );
                }
            }
        }
    }
    #[test]
    fn test_icc_of_tridiagonal_matrix_is_the_exact_cholesky_solve() {
        // No fill-in is dropped for a tridiagonal matrix, so M^-1 = A^-1.
        let n = 2000;
        let a = tridiagonal_csr(n, 2.5, -1.0);
        let icc = IccPreconditioner::new(&a);
        assert_eq!(icc.shift(), 0.0);
        let x_exact: Vec<f64> = (0..n).map(|i| ((i % 17) as f64 - 8.0) / 3.0).collect();
        let rhs = a.mul_vec(&x_exact);
        let z = icc.solve(&rhs);
        let worst = z
            .iter()
            .zip(x_exact.iter())
            .map(|(zi, xi)| (zi - xi).abs())
            .fold(0.0_f64, f64::max);
        assert!(worst < 1e-11, "max |z - x| = {worst}");
    }
    #[test]
    fn test_icc_of_tridiagonal_matrix_matches_dense_cholesky_entries() {
        let n = 6;
        let a = tridiagonal_csr(n, 2.0, -1.0);
        let icc = IccPreconditioner::new(&a);
        // Dense Cholesky of tridiag(-1, 2, -1): L_ii = sqrt((i + 2) / (i + 1)),
        // L_{i+1,i} = -1 / L_ii.
        for i in 0..n {
            let l_ii = ((i as f64 + 2.0) / (i as f64 + 1.0)).sqrt();
            let diag_pos = icc.row_ptr[i + 1] - 1;
            assert_eq!(icc.col_indices[diag_pos], i);
            assert!((icc.l_values[diag_pos] - l_ii).abs() < 1e-13, "L[{i},{i}]");
            if i + 1 < n {
                let next_row_start = icc.row_ptr[i + 1];
                assert_eq!(icc.col_indices[next_row_start], i);
                assert!((icc.l_values[next_row_start] + 1.0 / l_ii).abs() < 1e-13);
            }
        }
    }
    #[test]
    fn test_icc_factorization_of_large_poisson_matrix_is_linear_time() {
        // The column-by-all-rows search of the previous implementation was
        // O(n^2 * row_nnz): minutes to hours for n = 1e5.  The row-oriented
        // factorization takes milliseconds.  An optimised build must stay
        // within 100 ms; unoptimised test builds run the same loops an order
        // of magnitude slower on a machine shared with other test processes, so
        // there the bound only guards against a return to quadratic behaviour.
        let budget = if cfg!(debug_assertions) {
            std::time::Duration::from_secs(2)
        } else {
            std::time::Duration::from_millis(100)
        };
        let n = 100_000;
        let a = tridiagonal_csr(n, 2.0, -1.0);
        let started = std::time::Instant::now();
        let icc = IccPreconditioner::new(&a);
        let elapsed = started.elapsed();
        assert!(
            elapsed < budget,
            "IC(0) of a 1e5 tridiagonal matrix took {elapsed:?} (budget {budget:?})"
        );
        assert_eq!(icc.nnz(), 2 * n - 1);
        assert_eq!(icc.shift(), 0.0);
        // L_ii^2 = (i + 2) / (i + 1) exactly for the 1-D Poisson matrix.
        for i in [0, 1, 17, 4_999, n - 1] {
            let d = icc.l_values[icc.row_ptr[i + 1] - 1];
            let expected = ((i as f64 + 2.0) / (i as f64 + 1.0)).sqrt();
            assert!((d - expected).abs() < 1e-9, "L[{i},{i}] = {d}");
        }
    }
    /// The 4-cycle matrix with edge weights `t` and one negated edge: symmetric
    /// positive definite for `t < 1/sqrt(2)` (eigenvalues `1 +- sqrt(2) t`), but
    /// IC(0) drops the fill-in at (2, 1) and the last pivot becomes
    /// `1 - 2 t^2 / (1 - t^2) < 0` once `t^2 > 1/3`.
    fn frustrated_four_cycle(t: f64) -> CsrMatrix {
        CsrMatrix::from_triplets(
            4,
            4,
            &[
                (0, 0, 1.0),
                (1, 1, 1.0),
                (2, 2, 1.0),
                (3, 3, 1.0),
                (0, 1, t),
                (1, 0, t),
                (0, 2, t),
                (2, 0, t),
                (1, 3, t),
                (3, 1, t),
                (2, 3, -t),
                (3, 2, -t),
            ],
        )
    }
    #[test]
    fn test_icc_breakdown_is_cured_by_manteuffel_shift() {
        let t = 0.65;
        let a = frustrated_four_cycle(t);
        // The matrix is positive definite: its smallest eigenvalue is 1 - sqrt(2) t.
        let dense: Vec<Vec<f64>> = (0..4)
            .map(|i| (0..4).map(|j| a.get(i, j)).collect())
            .collect();
        let (eigenvalues, _) = crate::eigen::jacobi_eigen_dense(&dense, 4);
        let lambda_min = eigenvalues.iter().cloned().fold(f64::INFINITY, f64::min);
        assert!(
            (lambda_min - (1.0 - 2.0_f64.sqrt() * t)).abs() < 1e-12,
            "lambda_min = {lambda_min}"
        );
        // Plain IC(0) has the pivot 1 - 2 t^2 / (1 - t^2) = -0.463 < 0, so the
        // factorization must have been restarted with a shift.  The shifts
        // 1e-3 * 2^k fail for k <= 6 and the first success is 0.128.
        let icc = IccPreconditioner::new(&a);
        assert!(!icc.used_diagonal_fallback());
        assert!(
            (icc.shift() - 0.128).abs() < 1e-12,
            "shift = {}",
            icc.shift()
        );
        // The factor is a bounded Cholesky factor of A + shift * diag(A) on the
        // pattern, not the 1e-15 clamp that used to blow up the solve.
        let max_entry = icc.l_values.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
        let min_diag = (0..4)
            .map(|i| icc.l_values[icc.row_ptr[i + 1] - 1])
            .fold(f64::INFINITY, f64::min);
        assert!(max_entry < 2.0, "max |L| = {max_entry}");
        assert!(min_diag > 0.05, "min L_ii = {min_diag}");
        for i in 0..4 {
            for p in a.row_ptr[i]..a.row_ptr[i + 1] {
                let j = a.col_indices[p];
                if j > i {
                    continue;
                }
                let target = if i == j { 1.128 } else { a.values[p] };
                let reproduced = icc_product_entry(&icc, i, j);
                assert!(
                    (reproduced - target).abs() < 1e-12,
                    "[{i},{j}] {reproduced}"
                );
            }
        }
        // As a preconditioner it is symmetric positive definite and PCG
        // converges on the exact-arithmetic-terminating 4 x 4 system.
        let b = vec![1.0, -2.0, 0.5, 3.0];
        let (x, stats) = crate::solvers::pcg_icc(&a, &b, 50, 1e-10);
        assert!(stats.converged, "iterations = {}", stats.iterations);
        assert!(stats.iterations <= 8, "iterations = {}", stats.iterations);
        let ax = a.mul_vec(&x);
        for (axi, bi) in ax.iter().zip(b.iter()) {
            assert!((axi - bi).abs() < 1e-9);
        }
        let z = icc.solve(&b);
        let zb: f64 = z.iter().zip(b.iter()).map(|(zi, bi)| zi * bi).sum();
        assert!(zb > 0.0, "preconditioner must be positive definite: {zb}");
    }
    #[test]
    fn test_icc_indefinite_matrix_falls_back_to_diagonal_scaling() {
        // [[4, 10], [10, 9]] has determinant -64 < 0 (indefinite).  With the
        // diagonal scaled by 1 + s the second pivot is 9 (1 + s) - 25 / (1 + s),
        // positive only for s > 2/3, beyond the largest shift 0.512 tried.
        let a = CsrMatrix::from_triplets(
            2,
            2,
            &[(0, 0, 4.0), (0, 1, 10.0), (1, 0, 10.0), (1, 1, 9.0)],
        );
        let icc = IccPreconditioner::new(&a);
        assert!(icc.used_diagonal_fallback());
        assert!(icc.shift().is_infinite());
        let z = icc.solve(&[8.0, 27.0]);
        assert!((z[0] - 2.0).abs() < 1e-14, "z0 = {}", z[0]);
        assert!((z[1] - 3.0).abs() < 1e-14, "z1 = {}", z[1]);
    }
    #[test]
    fn test_icc_missing_or_non_positive_diagonal_does_not_poison_the_factor() {
        // Row 1 stores no diagonal entry at all.
        let a = CsrMatrix::from_triplets(3, 3, &[(0, 0, 4.0), (1, 0, 1.0), (2, 2, 9.0)]);
        let icc = IccPreconditioner::new(&a);
        assert!(icc.used_diagonal_fallback());
        let z = icc.solve(&[4.0, 5.0, 18.0]);
        assert!(z.iter().all(|v| v.is_finite()), "z = {z:?}");
        assert!((z[0] - 1.0).abs() < 1e-14);
        assert!(
            (z[1] - 5.0).abs() < 1e-14,
            "unit scaling for the empty diagonal"
        );
        assert!((z[2] - 2.0).abs() < 1e-14);
    }
    #[test]
    fn test_icc_accepts_unsorted_and_duplicate_csr_entries() {
        let sorted = tridiagonal_csr(4, 3.0, -1.0);
        // The same matrix with scrambled column order in every row, the diagonal
        // split into two summed halves, and the upper triangle stored as well.
        let scrambled = CsrMatrix {
            row_ptr: vec![0, 3, 7, 11, 14],
            col_indices: vec![1, 0, 0, 2, 1, 0, 1, 3, 2, 1, 2, 3, 2, 3],
            values: vec![
                -1.0, 1.5, 1.5, -1.0, 1.5, -1.0, 1.5, -1.0, 1.5, -1.0, 1.5, 1.5, -1.0, 1.5,
            ],
            nrows: 4,
            ncols: 4,
        };
        let reference = IccPreconditioner::new(&sorted);
        let icc = IccPreconditioner::new(&scrambled);
        assert_eq!(icc.row_ptr, reference.row_ptr);
        assert_eq!(icc.col_indices, reference.col_indices);
        for (got, want) in icc.l_values.iter().zip(reference.l_values.iter()) {
            assert!((got - want).abs() < 1e-14, "{got} vs {want}");
        }
    }
    #[test]
    fn test_icc_pcg_needs_far_fewer_iterations_than_plain_cg() {
        let a = poisson_2d_csr(24);
        let n = a.nrows;
        let b: Vec<f64> = (0..n).map(|i| 1.0 + (i % 5) as f64).collect();
        let (x, icc_stats) = crate::solvers::pcg_icc(&a, &b, 400, 1e-10);
        assert!(icc_stats.converged);
        let ax = a.mul_vec(&x);
        let resid = ax
            .iter()
            .zip(b.iter())
            .map(|(p, q)| (p - q).powi(2))
            .sum::<f64>()
            .sqrt();
        let b_norm = b.iter().map(|v| v * v).sum::<f64>().sqrt();
        assert!(
            resid / b_norm < 1e-9,
            "relative residual {}",
            resid / b_norm
        );
        let identity = vec![1.0; n];
        let (_, cg_stats) = crate::solvers::pcg_solve(&a, &b, &identity, 400, 1e-10);
        assert!(cg_stats.converged);
        // IC(0) on the 5-point Laplacian typically needs ~55-60 % of the CG
        // iterations (the condition number is reduced by a constant factor).
        assert!(
            10 * icc_stats.iterations < 7 * cg_stats.iterations,
            "ICC-PCG: {} iterations, CG: {}",
            icc_stats.iterations,
            cg_stats.iterations
        );
    }
    #[test]
    fn test_icc_pcg_on_assembled_elasticity_matrix_beats_jacobi_pcg() {
        use crate::assembly::{apply_dirichlet_bc, assemble_stiffness};
        use crate::constitutive::LinearElasticMaterial;
        use crate::mesh::TetrahedralMesh;
        // Unit cube of linear tetrahedra, clamped at z = 0 and sheared and
        // pressed on z = 1.  A tetrahedral elasticity matrix is not an
        // M-matrix, so this is a realistic case for IC(0).
        let mesh = TetrahedralMesh::generate_beam(1.0, 1.0, 1.0, 5, 5, 5);
        let material = LinearElasticMaterial::new(1.0e3, 0.3);
        let mut k = assemble_stiffness(&mesh, &material);
        let n = k.nrows;
        let mut f = vec![0.0_f64; n];
        let mut fixed: Vec<usize> = Vec::new();
        for (node, p) in mesh.nodes.iter().enumerate() {
            if p.z < 1.0e-9 {
                fixed.extend([3 * node, 3 * node + 1, 3 * node + 2]);
            } else if p.z > 1.0 - 1.0e-9 {
                f[3 * node] = 1.0;
                f[3 * node + 2] = -0.5;
            }
        }
        apply_dirichlet_bc(&mut k, &mut f, &fixed).expect("fixed DOFs are in range");
        let f_norm = f.iter().map(|v| v * v).sum::<f64>().sqrt();
        assert!(f_norm > 1.0, "the load must not vanish: {f_norm}");
        let icc = IccPreconditioner::new(&k);
        assert!(!icc.used_diagonal_fallback(), "shift = {}", icc.shift());
        assert!(icc.shift() < 0.1, "shift = {}", icc.shift());
        let (x, icc_stats) = crate::solvers::pcg_icc(&k, &f, 1000, 1e-10);
        assert!(icc_stats.converged, "ICC-PCG: {icc_stats:?}");
        let residual = k
            .mul_vec(&x)
            .iter()
            .zip(f.iter())
            .map(|(ax, b)| (ax - b).powi(2))
            .sum::<f64>()
            .sqrt();
        assert!(
            residual < 1e-8 * f_norm,
            "relative residual {}",
            residual / f_norm
        );
        let jacobi: Vec<f64> = (0..n).map(|i| 1.0 / k.get(i, i)).collect();
        let (_, jacobi_stats) = crate::solvers::pcg_solve(&k, &f, &jacobi, 3000, 1e-10);
        assert!(jacobi_stats.converged, "Jacobi-PCG: {jacobi_stats:?}");
        assert!(
            5 * icc_stats.iterations < 4 * jacobi_stats.iterations,
            "ICC-PCG: {} iterations, Jacobi-PCG: {}",
            icc_stats.iterations,
            jacobi_stats.iterations
        );
    }
    #[test]
    fn test_block_csr_identity_mul_vec() {
        let id3: [f64; 9] = [1., 0., 0., 0., 1., 0., 0., 0., 1.];
        let triplets = vec![(0usize, 0usize, id3), (1, 1, id3)];
        let bm = BlockCsrMatrix3::from_block_triplets(2, 2, &triplets);
        let x = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let y = bm.mul_vec(&x);
        for i in 0..6 {
            assert!(
                (y[i] - x[i]).abs() < 1e-12,
                "y[{i}] = {}, expected {}",
                y[i],
                x[i]
            );
        }
    }
    #[test]
    fn test_block_csr_sum_duplicate_blocks() {
        let blk1: [f64; 9] = [1., 0., 0., 0., 1., 0., 0., 0., 1.];
        let blk2: [f64; 9] = [2., 0., 0., 0., 2., 0., 0., 0., 2.];
        let triplets = vec![(0, 0, blk1), (0, 0, blk2)];
        let bm = BlockCsrMatrix3::from_block_triplets(1, 1, &triplets);
        assert_eq!(bm.n_blocks(), 1);
        let result = bm.get_block(0, 0);
        assert!(
            (result[0] - 3.0).abs() < 1e-12,
            "block[0,0] = {}",
            result[0]
        );
        assert!(
            (result[4] - 3.0).abs() < 1e-12,
            "block[1,1] = {}",
            result[4]
        );
    }
    #[test]
    fn test_block_csr_to_scalar_csr() {
        let blk: [f64; 9] = [1., 2., 3., 4., 5., 6., 7., 8., 9.];
        let bm = BlockCsrMatrix3::from_block_triplets(1, 1, &[(0, 0, blk)]);
        let scalar = bm.to_scalar_csr();
        assert_eq!(scalar.nrows, 3);
        assert_eq!(scalar.ncols, 3);
        for i in 0..3 {
            for j in 0..3 {
                let expected = blk[i * 3 + j];
                let got = scalar.get(i, j);
                assert!(
                    (got - expected).abs() < 1e-12,
                    "scalar[{i},{j}] = {got}, expected {expected}"
                );
            }
        }
    }
    #[test]
    fn test_block_csr_frobenius_norm() {
        let id3: [f64; 9] = [1., 0., 0., 0., 1., 0., 0., 0., 1.];
        let bm = BlockCsrMatrix3::from_block_triplets(1, 1, &[(0, 0, id3)]);
        assert!((bm.frobenius_norm() - 3.0_f64.sqrt()).abs() < 1e-12);
    }
    #[test]
    fn test_csr_norm_1() {
        let m =
            CsrMatrix::from_triplets(2, 2, &[(0, 0, 1.0), (0, 1, 2.0), (1, 0, 3.0), (1, 1, 4.0)]);
        assert!(
            (csr_norm_1(&m) - 6.0).abs() < 1e-12,
            "1-norm = {}",
            csr_norm_1(&m)
        );
    }
    #[test]
    fn test_csr_norm_inf() {
        let m =
            CsrMatrix::from_triplets(2, 2, &[(0, 0, 1.0), (0, 1, 2.0), (1, 0, 3.0), (1, 1, 4.0)]);
        assert!(
            (csr_norm_inf(&m) - 7.0).abs() < 1e-12,
            "inf-norm = {}",
            csr_norm_inf(&m)
        );
    }
    #[test]
    fn test_rcm_is_valid_permutation() {
        let triplets = vec![
            (0, 0, 4.0),
            (0, 1, -1.0),
            (0, 3, -1.0),
            (1, 0, -1.0),
            (1, 1, 4.0),
            (1, 2, -1.0),
            (2, 1, -1.0),
            (2, 2, 4.0),
            (2, 3, -1.0),
            (3, 0, -1.0),
            (3, 2, -1.0),
            (3, 3, 4.0),
        ];
        let m = CsrMatrix::from_triplets(4, 4, &triplets);
        let perm = reverse_cuthill_mckee(&m);
        assert_eq!(perm.len(), 4, "permutation must have n entries");
        let mut sorted = perm.clone();
        sorted.sort_unstable();
        assert_eq!(sorted, vec![0, 1, 2, 3], "permutation must be a bijection");
    }
    #[test]
    fn test_rcm_reduces_bandwidth() {
        let triplets = vec![
            (0, 0, 4.0),
            (0, 1, -1.0),
            (0, 2, -1.0),
            (0, 3, -1.0),
            (1, 0, -1.0),
            (1, 1, 4.0),
            (2, 0, -1.0),
            (2, 2, 4.0),
            (3, 0, -1.0),
            (3, 3, 4.0),
        ];
        let m = CsrMatrix::from_triplets(4, 4, &triplets);
        let id_perm: Vec<usize> = (0..4).collect();
        let bw_orig = bandwidth(&m, &id_perm);
        let perm = reverse_cuthill_mckee(&m);
        let bw_rcm = bandwidth(&m, &perm);
        assert!(
            bw_rcm <= bw_orig + 1,
            "RCM bandwidth ({bw_rcm}) should be <= original ({bw_orig}) + 1"
        );
    }
    #[test]
    fn test_to_upper_triangular() {
        let m = CsrMatrix::from_triplets(
            3,
            3,
            &[
                (0, 0, 1.0),
                (0, 1, 2.0),
                (0, 2, 3.0),
                (1, 0, 2.0),
                (1, 1, 4.0),
                (1, 2, 5.0),
                (2, 0, 3.0),
                (2, 1, 5.0),
                (2, 2, 6.0),
            ],
        );
        let upper = to_upper_triangular(&m);
        assert_eq!(upper.get(1, 0), 0.0, "lower entry should be 0");
        assert_eq!(upper.get(2, 0), 0.0, "lower entry should be 0");
        assert_eq!(upper.get(2, 1), 0.0, "lower entry should be 0");
        assert!((upper.get(0, 1) - 2.0).abs() < 1e-12);
        assert!((upper.get(0, 2) - 3.0).abs() < 1e-12);
        assert!((upper.get(1, 2) - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_from_upper_triangular_round_trip() {
        let upper = CsrMatrix::from_triplets(
            3,
            3,
            &[
                (0, 0, 1.0),
                (0, 1, 2.0),
                (0, 2, 3.0),
                (1, 1, 4.0),
                (1, 2, 5.0),
                (2, 2, 6.0),
            ],
        );
        let full = from_upper_triangular(&upper);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (full.get(i, j) - full.get(j, i)).abs() < 1e-12,
                    "not symmetric at ({i},{j})"
                );
            }
        }
        assert!((full.get(1, 0) - 2.0).abs() < 1e-12);
        assert!((full.get(2, 0) - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_is_symmetric_true() {
        let m =
            CsrMatrix::from_triplets(2, 2, &[(0, 0, 4.0), (0, 1, 1.0), (1, 0, 1.0), (1, 1, 3.0)]);
        assert!(is_symmetric(&m, 1e-12), "matrix should be symmetric");
    }
    #[test]
    fn test_is_symmetric_false() {
        let m =
            CsrMatrix::from_triplets(2, 2, &[(0, 0, 4.0), (0, 1, 2.0), (1, 0, 1.0), (1, 1, 3.0)]);
        assert!(!is_symmetric(&m, 1e-12), "matrix should not be symmetric");
    }
}
