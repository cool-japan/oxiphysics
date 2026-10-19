//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::collections::HashMap;

/// Incomplete Cholesky (ICC) preconditioner for symmetric positive definite
/// sparse matrices.
///
/// Computes a lower-triangular factor `L` such that `L L^T ≈ A`, maintaining
/// the sparsity pattern of the lower triangle of `A` (zero fill-in, "IC(0)").
/// `L L^T` reproduces `A` exactly on that pattern; for a tridiagonal `A` the
/// factor is the exact Cholesky factor.  This is the sparse counterpart of the
/// dense incomplete Cholesky; suitable for use as a preconditioner in PCG.
///
/// The factorization is row oriented ("up-looking"): row `i` is built from the
/// already finished rows `k < i` that it couples to, with every inner product a
/// sorted merge of two sparse rows, so the cost is
/// `O(sum_i sum_{k in row i} (nnz_i + nnz_k))`, linear in `nnz` for matrices
/// with a bounded row length.
///
/// # Breakdown handling
///
/// IC(0) can break down for a symmetric positive definite matrix that is not an
/// M- or H-matrix: dropping fill-in can make a pivot `A_ii - sum_k L_ik^2`
/// non-positive.  The factorization is then restarted on the diagonally shifted
/// matrix `A + α diag(A)` (Manteuffel, 1980) with `α = 1e-3, 2e-3, 4e-3, ...`
/// for at most ten shifts (`α <= 0.512`); the shift that finally succeeded is
/// reported by [`shift`](Self::shift).  If even the largest shift breaks down
/// the matrix is far from positive definite, and the preconditioner degrades to
/// the diagonal (Jacobi) scaling `L = diag(sqrt(A_ii))` (with `1` in place of any
/// non-positive diagonal entry) instead of returning an unusable factor; see
/// [`used_diagonal_fallback`](Self::used_diagonal_fallback).
#[derive(Debug, Clone)]
pub struct IccPreconditioner {
    /// Lower-triangular Cholesky factor stored in CSR format.  Every row holds
    /// its off-diagonal entries in ascending column order followed by the
    /// diagonal entry.
    pub(super) l_values: Vec<f64>,
    pub(super) row_ptr: Vec<usize>,
    pub(super) col_indices: Vec<usize>,
    pub(super) n: usize,
    /// Manteuffel shift `α` of the successful factorization.
    pub(super) shift: f64,
    /// `true` if the factor is the diagonal fallback.
    pub(super) diagonal_fallback: bool,
}
/// First Manteuffel shift tried after a plain IC(0) breakdown.
const ICC_INITIAL_SHIFT: f64 = 1.0e-3;
/// Number of shifted restarts (each doubling the shift) before the diagonal
/// fallback is used.
const ICC_MAX_SHIFT_TRIES: usize = 10;
/// A pivot not exceeding this fraction of its (shifted) diagonal entry is a
/// breakdown: the factor would otherwise contain entries that are larger by the
/// reciprocal of this fraction.
const ICC_PIVOT_TOLERANCE: f64 = 1.0e-10;
/// Inner product of two sparse vectors given by ascending index lists.
fn sorted_sparse_dot(idx_a: &[usize], val_a: &[f64], idx_b: &[usize], val_b: &[f64]) -> f64 {
    let mut ia = 0;
    let mut ib = 0;
    let mut sum = 0.0;
    while ia < idx_a.len() && ib < idx_b.len() {
        match idx_a[ia].cmp(&idx_b[ib]) {
            std::cmp::Ordering::Equal => {
                sum += val_a[ia] * val_b[ib];
                ia += 1;
                ib += 1;
            }
            std::cmp::Ordering::Less => ia += 1,
            std::cmp::Ordering::Greater => ib += 1,
        }
    }
    sum
}
impl IccPreconditioner {
    /// Compute the ICC factorization of a symmetric CSR matrix `a`.
    ///
    /// Only the lower-triangular pattern (including diagonal) of `a` is used.
    /// Fill-in is dropped (ICC(0) strategy).  Column indices within a row of `a`
    /// need not be sorted and duplicate entries are summed.  See the type-level
    /// documentation for the handling of pivot breakdown.
    ///
    /// # Panics
    ///
    /// Panics if `a` is not square.
    pub fn new(a: &CsrMatrix) -> Self {
        assert_eq!(a.nrows, a.ncols, "ICC requires a square matrix");
        let n = a.nrows;
        let (row_ptr, col_indices, a_lower) = Self::lower_triangle(a);
        let mut shift = 0.0_f64;
        for attempt in 0..=ICC_MAX_SHIFT_TRIES {
            if let Some(l_values) = Self::factorize(&row_ptr, &col_indices, &a_lower, shift) {
                return IccPreconditioner {
                    l_values,
                    row_ptr,
                    col_indices,
                    n,
                    shift,
                    diagonal_fallback: false,
                };
            }
            shift = if attempt == 0 {
                ICC_INITIAL_SHIFT
            } else {
                2.0 * shift
            };
        }
        let mut l_values = vec![0.0_f64; a_lower.len()];
        for i in 0..n {
            let diag_pos = row_ptr[i + 1] - 1;
            let d = a_lower[diag_pos];
            l_values[diag_pos] = if d.is_finite() && d > 0.0 {
                d.sqrt()
            } else {
                1.0
            };
        }
        IccPreconditioner {
            l_values,
            row_ptr,
            col_indices,
            n,
            shift: f64::INFINITY,
            diagonal_fallback: true,
        }
    }
    /// Diagonal shift `α` of the factorization: the factor satisfies
    /// `L L^T ≈ A + α diag(A)`.
    ///
    /// `0.0` for a plain IC(0) factor; a positive Manteuffel shift when IC(0)
    /// broke down; `f64::INFINITY` for the diagonal fallback.
    pub fn shift(&self) -> f64 {
        self.shift
    }
    /// `true` if the shifted factorizations all broke down and the
    /// preconditioner is the diagonal (Jacobi) scaling.
    pub fn used_diagonal_fallback(&self) -> bool {
        self.diagonal_fallback
    }
    /// Lower triangle of `a` (diagonal included) as CSR arrays
    /// `(row_ptr, col_indices, values)`.  Each row is sorted by column with
    /// duplicates summed and the diagonal entry present (as an explicit zero if
    /// `a` does not store it), so the diagonal is the last entry of every row.
    fn lower_triangle(a: &CsrMatrix) -> (Vec<usize>, Vec<usize>, Vec<f64>) {
        let n = a.nrows;
        let mut row_ptr = Vec::with_capacity(n + 1);
        let mut col_indices: Vec<usize> = Vec::with_capacity(a.nnz() / 2 + n);
        let mut values: Vec<f64> = Vec::with_capacity(a.nnz() / 2 + n);
        row_ptr.push(0);
        let mut entries: Vec<(usize, f64)> = Vec::new();
        for row in 0..n {
            entries.clear();
            for idx in a.row_ptr[row]..a.row_ptr[row + 1] {
                let col = a.col_indices[idx];
                if col <= row {
                    entries.push((col, a.values[idx]));
                }
            }
            entries.sort_by_key(|&(col, _)| col);
            let row_start = col_indices.len();
            for &(col, value) in &entries {
                if col_indices.len() > row_start && col_indices.last() == Some(&col) {
                    if let Some(last) = values.last_mut() {
                        *last += value;
                    }
                } else {
                    col_indices.push(col);
                    values.push(value);
                }
            }
            if col_indices.last() != Some(&row) || col_indices.len() == row_start {
                col_indices.push(row);
                values.push(0.0);
            }
            row_ptr.push(col_indices.len());
        }
        (row_ptr, col_indices, values)
    }
    /// Row-oriented IC(0) of `A + shift * diag(A)` on the pattern given by the
    /// lower-triangle CSR arrays.  Returns `None` on pivot breakdown.
    fn factorize(
        row_ptr: &[usize],
        col_indices: &[usize],
        a_lower: &[f64],
        shift: f64,
    ) -> Option<Vec<f64>> {
        let n = row_ptr.len().saturating_sub(1);
        let mut l = vec![0.0_f64; a_lower.len()];
        for i in 0..n {
            let row_start = row_ptr[i];
            let diag_pos = row_ptr[i + 1] - 1;
            for p in row_start..diag_pos {
                let k = col_indices[p];
                let k_start = row_ptr[k];
                let k_diag = row_ptr[k + 1] - 1;
                let dot = sorted_sparse_dot(
                    &col_indices[row_start..p],
                    &l[row_start..p],
                    &col_indices[k_start..k_diag],
                    &l[k_start..k_diag],
                );
                l[p] = (a_lower[p] - dot) / l[k_diag];
            }
            let sum_sq: f64 = l[row_start..diag_pos].iter().map(|&v| v * v).sum();
            let a_ii = (1.0 + shift) * a_lower[diag_pos];
            let pivot = a_ii - sum_sq;
            if !(a_ii > 0.0 && pivot > ICC_PIVOT_TOLERANCE * a_ii) {
                return None;
            }
            l[diag_pos] = pivot.sqrt();
        }
        Some(l)
    }
    /// Apply the ICC preconditioner: solve `(L L^T) z = r`.
    ///
    /// Performs two triangular solves: `L y = r` (forward), then `L^T z = y`
    /// (backward).
    pub fn solve(&self, rhs: &[f64]) -> Vec<f64> {
        let n = self.n;
        assert_eq!(rhs.len(), n);
        let mut y = rhs.to_vec();
        for i in 0..n {
            let start = self.row_ptr[i];
            let end = self.row_ptr[i + 1];
            let mut diag = 1.0;
            for p in start..end {
                let j = self.col_indices[p];
                if j == i {
                    diag = self.l_values[p];
                } else if j < i {
                    y[i] -= self.l_values[p] * y[j];
                }
            }
            if diag.abs() > 1e-60 {
                y[i] /= diag;
            }
        }
        let mut z = y.clone();
        for i in (0..n).rev() {
            let start = self.row_ptr[i];
            let end = self.row_ptr[i + 1];
            let mut diag = 1.0;
            for p in start..end {
                let j = self.col_indices[p];
                if j == i {
                    diag = self.l_values[p];
                }
            }
            if diag.abs() > 1e-60 {
                z[i] /= diag;
            }
            for p in start..end {
                let j = self.col_indices[p];
                if j < i {
                    z[j] -= self.l_values[p] * z[i];
                }
            }
        }
        z
    }
    /// Return the number of non-zero entries in the factor L.
    pub fn nnz(&self) -> usize {
        self.l_values.len()
    }
}
/// A node in a 2-D quadtree used for adaptive mesh refinement.
///
/// Each node represents a rectangular cell in the mesh.  Leaf nodes
/// correspond to actual mesh cells; internal nodes have been refined.
#[derive(Debug, Clone)]
pub struct QuadTreeNode {
    /// Lower-left x coordinate.
    pub x0: f64,
    /// Lower-left y coordinate.
    pub y0: f64,
    /// Cell width.
    pub width: f64,
    /// Cell height.
    pub height: f64,
    /// Refinement level (0 = root).
    pub level: u32,
    /// Per-element error indicator (set by the caller before refinement).
    pub error_indicator: f64,
    /// Children (SW, SE, NW, NE) if this node has been refined.
    pub children: Option<Box<[QuadTreeNode; 4]>>,
    /// Unique cell index (assigned during leaf enumeration).
    pub cell_id: usize,
}
impl QuadTreeNode {
    /// Create a root node covering `[x0, x0+width] × [y0, y0+height]`.
    pub fn new_root(x0: f64, y0: f64, width: f64, height: f64) -> Self {
        Self {
            x0,
            y0,
            width,
            height,
            level: 0,
            error_indicator: 0.0,
            children: None,
            cell_id: 0,
        }
    }
    /// Create a child node.
    fn new_child(x0: f64, y0: f64, width: f64, height: f64, level: u32) -> Self {
        Self {
            x0,
            y0,
            width,
            height,
            level,
            error_indicator: 0.0,
            children: None,
            cell_id: 0,
        }
    }
    /// Return the cell width (same as `self.width`).
    pub fn cell_width(&self) -> f64 {
        self.width
    }
    /// Return the cell height (same as `self.height`).
    pub fn cell_height(&self) -> f64 {
        self.height
    }
    /// Return `true` if this node has no children (is a leaf).
    pub fn is_leaf(&self) -> bool {
        self.children.is_none()
    }
    /// Refine this node into four child cells (SW, SE, NW, NE).
    ///
    /// Does nothing if the node is already refined.
    pub fn refine(&mut self) {
        if self.children.is_some() {
            return;
        }
        let hw = self.width * 0.5;
        let hh = self.height * 0.5;
        let lev = self.level + 1;
        let sw = Self::new_child(self.x0, self.y0, hw, hh, lev);
        let se = Self::new_child(self.x0 + hw, self.y0, hw, hh, lev);
        let nw = Self::new_child(self.x0, self.y0 + hh, hw, hh, lev);
        let ne = Self::new_child(self.x0 + hw, self.y0 + hh, hw, hh, lev);
        self.children = Some(Box::new([sw, se, nw, ne]));
    }
    /// Count the total number of leaf nodes in this subtree.
    pub fn leaf_count(&self) -> usize {
        match &self.children {
            None => 1,
            Some(ch) => ch.iter().map(|c| c.leaf_count()).sum(),
        }
    }
    /// Enumerate leaf nodes and assign sequential cell IDs.
    ///
    /// Returns the number of leaves assigned.
    pub fn enumerate_leaves(&mut self, start_id: usize) -> usize {
        match self.children.as_mut() {
            None => {
                self.cell_id = start_id;
                start_id + 1
            }
            Some(ch) => {
                let mut id = start_id;
                for child in ch.iter_mut() {
                    id = child.enumerate_leaves(id);
                }
                id
            }
        }
    }
    /// Collect all leaf nodes as immutable references.
    pub fn collect_leaves<'a>(&'a self, leaves: &mut Vec<&'a QuadTreeNode>) {
        match &self.children {
            None => leaves.push(self),
            Some(ch) => {
                for child in ch.iter() {
                    child.collect_leaves(leaves);
                }
            }
        }
    }
    /// Centroid x-coordinate.
    pub fn cx(&self) -> f64 {
        self.x0 + 0.5 * self.width
    }
    /// Centroid y-coordinate.
    pub fn cy(&self) -> f64 {
        self.y0 + 0.5 * self.height
    }
}
/// Incomplete LU factorization with level-of-fill `k` (ILU(k)).
///
/// Level k = 0 reproduces ILU(0) (no fill).  Level k = 1 allows one level
/// of fill beyond the original pattern, and so on.  Higher levels produce
/// better preconditioners at the cost of more memory.
///
/// The factorization is stored in CSR format with an extended sparsity pattern
/// computed via symbolic analysis.
#[derive(Debug, Clone)]
pub struct IlukPreconditioner {
    /// L and U factors in combined CSR (L strictly lower, U includes diagonal).
    pub(super) lu_values: Vec<f64>,
    pub(super) row_ptr: Vec<usize>,
    pub(super) col_indices: Vec<usize>,
    /// Level-of-fill for each non-zero.
    pub(super) fill_levels: Vec<u32>,
    pub(super) n: usize,
    /// Fill level parameter k.
    pub(super) k: u32,
}
impl IlukPreconditioner {
    /// Compute the ILU(k) factorization of matrix `a` with fill level `k`.
    ///
    /// For k=0 this matches ILU(0).  For k≥1 additional fill-in entries
    /// are allowed in the sparsity pattern.
    pub fn new(a: &CsrMatrix, k: u32) -> Self {
        assert_eq!(a.nrows, a.ncols, "matrix must be square");
        let n = a.nrows;
        let mut pattern: Vec<Vec<(usize, u32)>> = (0..n)
            .map(|row| {
                let start = a.row_ptr[row];
                let end = a.row_ptr[row + 1];
                let mut row_pat: Vec<(usize, u32)> =
                    (start..end).map(|idx| (a.col_indices[idx], 0)).collect();
                row_pat.sort_by_key(|&(c, _)| c);
                row_pat
            })
            .collect();
        for i in 1..n {
            let mut j = 0;
            while j < pattern[i].len() {
                let (col, lev_ij) = pattern[i][j];
                if col >= i {
                    break;
                }
                let row_p_copy: Vec<(usize, u32)> = pattern[col].clone();
                for &(q, lev_pq) in &row_p_copy {
                    if q <= col {
                        continue;
                    }
                    let new_lev = lev_ij + lev_pq + 1;
                    if new_lev <= k {
                        match pattern[i].binary_search_by_key(&q, |&(c, _)| c) {
                            Ok(pos) => {
                                if pattern[i][pos].1 > new_lev {
                                    pattern[i][pos].1 = new_lev;
                                }
                            }
                            Err(pos) => {
                                pattern[i].insert(pos, (q, new_lev));
                            }
                        }
                    }
                }
                j += 1;
            }
        }
        let mut row_ptr = vec![0usize; n + 1];
        for i in 0..n {
            row_ptr[i + 1] = row_ptr[i] + pattern[i].len();
        }
        let nnz = row_ptr[n];
        let mut col_indices = vec![0usize; nnz];
        let mut fill_levels = vec![0u32; nnz];
        let mut lu_values = vec![0.0f64; nnz];
        for i in 0..n {
            let start = row_ptr[i];
            for (k_idx, &(col, lev)) in pattern[i].iter().enumerate() {
                col_indices[start + k_idx] = col;
                fill_levels[start + k_idx] = lev;
            }
        }
        for row in 0..n {
            let a_start = a.row_ptr[row];
            let a_end = a.row_ptr[row + 1];
            for a_idx in a_start..a_end {
                let col = a.col_indices[a_idx];
                let lu_start = row_ptr[row];
                let lu_end = row_ptr[row + 1];
                if let Ok(offset) = col_indices[lu_start..lu_end].binary_search(&col) {
                    lu_values[lu_start + offset] = a.values[a_idx];
                }
            }
        }
        for i in 1..n {
            let row_start = row_ptr[i];
            let row_end = row_ptr[i + 1];
            let lower_cols: Vec<usize> = (row_start..row_end)
                .map(|p| col_indices[p])
                .take_while(|&c| c < i)
                .collect();
            for &kk in &lower_cols {
                let k_start = row_ptr[kk];
                let k_end = row_ptr[kk + 1];
                let diag_k = match col_indices[k_start..k_end].binary_search(&kk) {
                    Ok(off) => lu_values[k_start + off],
                    Err(_) => 0.0,
                };
                if diag_k.abs() < 1e-60 {
                    continue;
                }
                let p_ik = col_indices[row_start..row_end]
                    .binary_search(&kk)
                    .map(|off| row_start + off)
                    .unwrap_or(usize::MAX);
                if p_ik == usize::MAX {
                    continue;
                }
                lu_values[p_ik] /= diag_k;
                let factor = lu_values[p_ik];
                for k_idx in k_start..k_end {
                    let j = col_indices[k_idx];
                    if j <= kk {
                        continue;
                    }
                    if let Ok(off) = col_indices[row_start..row_end].binary_search(&j) {
                        lu_values[row_start + off] -= factor * lu_values[k_idx];
                    }
                }
            }
        }
        IlukPreconditioner {
            lu_values,
            row_ptr,
            col_indices,
            fill_levels,
            n,
            k,
        }
    }
    /// Apply the ILU(k) preconditioner: solve (LU) z = r.
    pub fn solve(&self, rhs: &[f64]) -> Vec<f64> {
        assert_eq!(rhs.len(), self.n);
        let n = self.n;
        let mut y = rhs.to_vec();
        for i in 0..n {
            let start = self.row_ptr[i];
            let end = self.row_ptr[i + 1];
            for p in start..end {
                let j = self.col_indices[p];
                if j >= i {
                    break;
                }
                y[i] -= self.lu_values[p] * y[j];
            }
        }
        for i in (0..n).rev() {
            let start = self.row_ptr[i];
            let end = self.row_ptr[i + 1];
            let mut diag = 1.0;
            for p in start..end {
                let j = self.col_indices[p];
                if j == i {
                    diag = self.lu_values[p];
                } else if j > i {
                    y[i] -= self.lu_values[p] * y[j];
                }
            }
            if diag.abs() > 1e-60 {
                y[i] /= diag;
            }
        }
        y
    }
    /// Return the fill level parameter used for this factorization.
    pub fn fill_level(&self) -> u32 {
        self.k
    }
    /// Return the number of non-zeros in the extended sparsity pattern.
    pub fn nnz(&self) -> usize {
        self.lu_values.len()
    }
    /// Return the fill-in count: entries added beyond the original pattern.
    pub fn fill_in_count(&self) -> usize {
        self.fill_levels.iter().filter(|&&lev| lev > 0).count()
    }
}
/// Block sparse matrix where each scalar entry is replaced by a 3×3 dense
/// block.  Suitable for 3-D elasticity problems where DOFs come in groups of
/// three per node.
///
/// The block structure is stored in CSR format on the block level.  Each
/// block `(i, j)` holds a flat 9-element row-major 3×3 sub-matrix.
#[derive(Debug, Clone)]
pub struct BlockCsrMatrix3 {
    /// Row pointer array (block-level, length = n_block_rows + 1).
    pub row_ptr: Vec<usize>,
    /// Block column indices (block-level).
    pub col_indices: Vec<usize>,
    /// Block values: each entry is a 9-element \[f64; 9\] (row-major 3×3 block).
    pub blocks: Vec<[f64; 9]>,
    /// Number of block rows.
    pub n_block_rows: usize,
    /// Number of block columns.
    pub n_block_cols: usize,
}
impl BlockCsrMatrix3 {
    /// Create an empty block CSR matrix.
    pub fn new(n_block_rows: usize, n_block_cols: usize) -> Self {
        Self {
            row_ptr: vec![0; n_block_rows + 1],
            col_indices: Vec::new(),
            blocks: Vec::new(),
            n_block_rows,
            n_block_cols,
        }
    }
    /// Build from a list of block-level triplets `(block_row, block_col, block_3x3)`.
    ///
    /// Duplicate block entries at the same `(block_row, block_col)` are summed
    /// element-wise.
    pub fn from_block_triplets(
        n_block_rows: usize,
        n_block_cols: usize,
        triplets: &[(usize, usize, [f64; 9])],
    ) -> Self {
        let mut map: std::collections::HashMap<(usize, usize), [f64; 9]> =
            std::collections::HashMap::new();
        for &(r, c, ref blk) in triplets {
            let entry = map.entry((r, c)).or_insert([0.0; 9]);
            for k in 0..9 {
                entry[k] += blk[k];
            }
        }
        let mut entries: Vec<((usize, usize), [f64; 9])> = map.into_iter().collect();
        entries.sort_by_key(|&((r, c), _)| (r, c));
        let mut row_ptr = vec![0usize; n_block_rows + 1];
        let mut col_indices = Vec::new();
        let mut blocks = Vec::new();
        for &((r, c), ref blk) in &entries {
            row_ptr[r + 1] += 1;
            col_indices.push(c);
            blocks.push(*blk);
        }
        for i in 1..=n_block_rows {
            row_ptr[i] += row_ptr[i - 1];
        }
        Self {
            row_ptr,
            col_indices,
            blocks,
            n_block_rows,
            n_block_cols,
        }
    }
    /// Multiply by a dense vector of length `3 * n_block_cols`.
    ///
    /// Returns a dense vector of length `3 * n_block_rows`.
    pub fn mul_vec(&self, x: &[f64]) -> Vec<f64> {
        assert_eq!(x.len(), self.n_block_cols * 3);
        let mut y = vec![0.0f64; self.n_block_rows * 3];
        for br in 0..self.n_block_rows {
            let start = self.row_ptr[br];
            let end = self.row_ptr[br + 1];
            for bidx in start..end {
                let bc = self.col_indices[bidx];
                let blk = &self.blocks[bidx];
                for i in 0..3 {
                    for j in 0..3 {
                        y[br * 3 + i] += blk[i * 3 + j] * x[bc * 3 + j];
                    }
                }
            }
        }
        y
    }
    /// Return the number of stored blocks.
    pub fn n_blocks(&self) -> usize {
        self.blocks.len()
    }
    /// Get a block at `(block_row, block_col)`.  Returns a zero block if not
    /// stored.
    pub fn get_block(&self, block_row: usize, block_col: usize) -> [f64; 9] {
        let start = self.row_ptr[block_row];
        let end = self.row_ptr[block_row + 1];
        for idx in start..end {
            if self.col_indices[idx] == block_col {
                return self.blocks[idx];
            }
        }
        [0.0; 9]
    }
    /// Frobenius norm of the entire block matrix.
    pub fn frobenius_norm(&self) -> f64 {
        self.blocks
            .iter()
            .flat_map(|b| b.iter())
            .map(|v| v * v)
            .sum::<f64>()
            .sqrt()
    }
    /// Convert to a scalar CSR matrix (expand each 3×3 block into scalar entries).
    pub fn to_scalar_csr(&self) -> CsrMatrix {
        let n_scalar_rows = self.n_block_rows * 3;
        let n_scalar_cols = self.n_block_cols * 3;
        let mut triplets = Vec::with_capacity(self.blocks.len() * 9);
        for br in 0..self.n_block_rows {
            let start = self.row_ptr[br];
            let end = self.row_ptr[br + 1];
            for bidx in start..end {
                let bc = self.col_indices[bidx];
                let blk = &self.blocks[bidx];
                for i in 0..3 {
                    for j in 0..3 {
                        let v = blk[i * 3 + j];
                        if v.abs() > 1e-30 {
                            triplets.push((br * 3 + i, bc * 3 + j, v));
                        }
                    }
                }
            }
        }
        CsrMatrix::from_triplets(n_scalar_rows, n_scalar_cols, &triplets)
    }
}
/// Compressed Sparse Row (CSR) matrix.
///
/// Stores a sparse matrix in CSR format with row pointers, column indices,
/// and values arrays. Efficient for matrix-vector multiplication and row access.
#[derive(Debug, Clone)]
pub struct CsrMatrix {
    /// Row pointer array (length = nrows + 1).
    pub row_ptr: Vec<usize>,
    /// Column indices for non-zero entries.
    pub col_indices: Vec<usize>,
    /// Non-zero values.
    pub values: Vec<f64>,
    /// Number of rows.
    pub nrows: usize,
    /// Number of columns.
    pub ncols: usize,
}
impl CsrMatrix {
    /// Create a new empty CSR matrix with the given dimensions.
    pub fn new(nrows: usize, ncols: usize) -> Self {
        Self {
            row_ptr: vec![0; nrows + 1],
            col_indices: Vec::new(),
            values: Vec::new(),
            nrows,
            ncols,
        }
    }
    /// Build a CSR matrix from coordinate (triplet) format.
    ///
    /// Duplicate entries at the same `(row, col)` are summed together.
    ///
    /// # Panics
    ///
    /// Panics if any index is out of bounds.
    pub fn from_triplets(nrows: usize, ncols: usize, triplets: &[(usize, usize, f64)]) -> Self {
        let mut map: HashMap<(usize, usize), f64> = HashMap::new();
        for &(r, c, v) in triplets {
            assert!(r < nrows, "row index {r} out of bounds for {nrows} rows");
            assert!(c < ncols, "col index {c} out of bounds for {ncols} cols");
            *map.entry((r, c)).or_insert(0.0) += v;
        }
        let mut entries: Vec<((usize, usize), f64)> = map.into_iter().collect();
        entries.sort_by_key(|&((r, c), _)| (r, c));
        let mut row_ptr = vec![0usize; nrows + 1];
        let mut col_indices = Vec::with_capacity(entries.len());
        let mut values = Vec::with_capacity(entries.len());
        for &((r, c), v) in &entries {
            row_ptr[r + 1] += 1;
            col_indices.push(c);
            values.push(v);
        }
        for i in 1..=nrows {
            row_ptr[i] += row_ptr[i - 1];
        }
        Self {
            row_ptr,
            col_indices,
            values,
            nrows,
            ncols,
        }
    }
    /// Get the value at `(row, col)`. Returns 0.0 if the entry is not stored.
    ///
    /// # Panics
    ///
    /// Panics if indices are out of bounds.
    pub fn get(&self, row: usize, col: usize) -> f64 {
        assert!(row < self.nrows);
        assert!(col < self.ncols);
        let start = self.row_ptr[row];
        let end = self.row_ptr[row + 1];
        for idx in start..end {
            if self.col_indices[idx] == col {
                return self.values[idx];
            }
        }
        0.0
    }
    /// Set the value at `(row, col)`. If the entry exists, update it;
    /// otherwise insert a new entry.
    ///
    /// **Note:** Inserting new entries is O(nnz) in the worst case because the
    /// arrays must be shifted. Prefer [`from_triplets`](Self::from_triplets) for
    /// bulk construction.
    ///
    /// # Panics
    ///
    /// Panics if indices are out of bounds.
    pub fn set(&mut self, row: usize, col: usize, value: f64) {
        assert!(row < self.nrows);
        assert!(col < self.ncols);
        let start = self.row_ptr[row];
        let end = self.row_ptr[row + 1];
        for idx in start..end {
            if self.col_indices[idx] == col {
                self.values[idx] = value;
                return;
            }
        }
        let mut insert_pos = start;
        while insert_pos < end && self.col_indices[insert_pos] < col {
            insert_pos += 1;
        }
        self.col_indices.insert(insert_pos, col);
        self.values.insert(insert_pos, value);
        for r in (row + 1)..=self.nrows {
            self.row_ptr[r] += 1;
        }
    }
    /// Add `value` to the entry at `(row, col)`. If the entry does not exist,
    /// it is created with the given value. This is the primary method used
    /// during stiffness matrix assembly.
    ///
    /// # Panics
    ///
    /// Panics if indices are out of bounds.
    pub fn add_to(&mut self, row: usize, col: usize, value: f64) {
        assert!(row < self.nrows);
        assert!(col < self.ncols);
        let start = self.row_ptr[row];
        let end = self.row_ptr[row + 1];
        for idx in start..end {
            if self.col_indices[idx] == col {
                self.values[idx] += value;
                return;
            }
        }
        let mut insert_pos = start;
        while insert_pos < end && self.col_indices[insert_pos] < col {
            insert_pos += 1;
        }
        self.col_indices.insert(insert_pos, col);
        self.values.insert(insert_pos, value);
        for r in (row + 1)..=self.nrows {
            self.row_ptr[r] += 1;
        }
    }
    /// Multiply this matrix by a dense vector: `y = A * x`.
    ///
    /// # Panics
    ///
    /// Panics if `x.len() != self.ncols`.
    pub fn mul_vec(&self, x: &[f64]) -> Vec<f64> {
        assert_eq!(x.len(), self.ncols, "vector length must equal ncols");
        let mut y = vec![0.0; self.nrows];
        for (row, y_row) in y.iter_mut().enumerate().take(self.nrows) {
            let start = self.row_ptr[row];
            let end = self.row_ptr[row + 1];
            let mut sum = 0.0;
            for idx in start..end {
                sum += self.values[idx] * x[self.col_indices[idx]];
            }
            *y_row = sum;
        }
        y
    }
    /// Return the number of stored non-zero entries.
    pub fn nnz(&self) -> usize {
        self.values.len()
    }
    /// Get the diagonal element at `(i, i)`.
    pub fn diagonal(&self, i: usize) -> f64 {
        self.get(i, i)
    }
    /// Transpose this CSR matrix, returning a new CSR matrix.
    pub fn transpose(&self) -> CsrMatrix {
        let mut triplets = Vec::with_capacity(self.nnz());
        for row in 0..self.nrows {
            let start = self.row_ptr[row];
            let end = self.row_ptr[row + 1];
            for idx in start..end {
                triplets.push((self.col_indices[idx], row, self.values[idx]));
            }
        }
        CsrMatrix::from_triplets(self.ncols, self.nrows, &triplets)
    }
    /// Add two CSR matrices: C = A + B.
    ///
    /// Both matrices must have the same dimensions.
    pub fn add(&self, other: &CsrMatrix) -> CsrMatrix {
        assert_eq!(self.nrows, other.nrows, "row dimensions must match");
        assert_eq!(self.ncols, other.ncols, "col dimensions must match");
        let mut triplets = Vec::with_capacity(self.nnz() + other.nnz());
        for row in 0..self.nrows {
            let start = self.row_ptr[row];
            let end = self.row_ptr[row + 1];
            for idx in start..end {
                triplets.push((row, self.col_indices[idx], self.values[idx]));
            }
        }
        for row in 0..other.nrows {
            let start = other.row_ptr[row];
            let end = other.row_ptr[row + 1];
            for idx in start..end {
                triplets.push((row, other.col_indices[idx], other.values[idx]));
            }
        }
        CsrMatrix::from_triplets(self.nrows, self.ncols, &triplets)
    }
    /// Scale all values by a scalar: A *= alpha.
    pub fn scale(&mut self, alpha: f64) {
        for v in self.values.iter_mut() {
            *v *= alpha;
        }
    }
    /// Return a scaled copy: B = alpha * A.
    pub fn scaled(&self, alpha: f64) -> CsrMatrix {
        let mut result = self.clone();
        result.scale(alpha);
        result
    }
    /// Optimized sparse matrix-vector multiply using row-based access.
    ///
    /// Same as `mul_vec` but with explicit prefetch-friendly ordering.
    /// y = alpha * A * x + beta * y
    pub fn mul_vec_axpby(&self, x: &[f64], y: &mut [f64], alpha: f64, beta: f64) {
        assert_eq!(x.len(), self.ncols);
        assert_eq!(y.len(), self.nrows);
        for (row, y_row) in y.iter_mut().enumerate() {
            let start = self.row_ptr[row];
            let end = self.row_ptr[row + 1];
            let mut sum = 0.0;
            for idx in start..end {
                sum += self.values[idx] * x[self.col_indices[idx]];
            }
            *y_row = alpha * sum + beta * *y_row;
        }
    }
    /// Symmetric sparse matrix-vector multiply (upper triangle only).
    ///
    /// Assumes only the upper triangle is stored. Computes y = A * x
    /// using symmetry: a_ij contributes to both y\[i\] and y\[j\].
    pub fn symmetric_mul_vec(&self, x: &[f64]) -> Vec<f64> {
        assert_eq!(x.len(), self.ncols);
        assert_eq!(
            self.nrows, self.ncols,
            "matrix must be square for symmetric multiply"
        );
        let n = self.nrows;
        let mut y = vec![0.0; n];
        for row in 0..n {
            let start = self.row_ptr[row];
            let end = self.row_ptr[row + 1];
            for idx in start..end {
                let col = self.col_indices[idx];
                let val = self.values[idx];
                y[row] += val * x[col];
                if col != row {
                    y[col] += val * x[row];
                }
            }
        }
        y
    }
    /// Extract diagonal as a vector.
    pub fn diagonal_vec(&self) -> Vec<f64> {
        let n = self.nrows.min(self.ncols);
        (0..n).map(|i| self.get(i, i)).collect()
    }
    /// Convert this CSR matrix to CSC format.
    pub fn to_csc(&self) -> CscMatrix {
        let mut triplets = Vec::with_capacity(self.nnz());
        for row in 0..self.nrows {
            let start = self.row_ptr[row];
            let end = self.row_ptr[row + 1];
            for idx in start..end {
                triplets.push((row, self.col_indices[idx], self.values[idx]));
            }
        }
        CscMatrix::from_triplets(self.nrows, self.ncols, &triplets)
    }
    /// Frobenius norm: ||A||_F = sqrt(sum of a_ij^2).
    pub fn frobenius_norm(&self) -> f64 {
        self.values.iter().map(|v| v * v).sum::<f64>().sqrt()
    }
}
/// Thin wrapper around `Vec`f64` for sparse/dense vector operations.
#[derive(Debug, Clone)]
pub struct SparseVector {
    /// The underlying dense data.
    pub data: Vec<f64>,
}
impl SparseVector {
    /// Create a new zero vector of the given length.
    pub fn new(len: usize) -> Self {
        Self {
            data: vec![0.0; len],
        }
    }
    /// Create a sparse vector from existing data.
    pub fn from_vec(data: Vec<f64>) -> Self {
        Self { data }
    }
    /// Return the length of the vector.
    pub fn len(&self) -> usize {
        self.data.len()
    }
    /// Check if the vector is empty.
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }
    /// Compute the dot product with another vector.
    pub fn dot(&self, other: &Self) -> f64 {
        self.data
            .iter()
            .zip(other.data.iter())
            .map(|(a, b)| a * b)
            .sum()
    }
    /// Compute the Euclidean norm.
    pub fn norm(&self) -> f64 {
        self.data.iter().map(|x| x * x).sum::<f64>().sqrt()
    }
    /// AXPY: self = alpha * self + beta * other
    pub fn axpby(&mut self, alpha: f64, other: &SparseVector, beta: f64) {
        assert_eq!(self.data.len(), other.data.len());
        for (a, b) in self.data.iter_mut().zip(other.data.iter()) {
            *a = alpha * *a + beta * *b;
        }
    }
    /// Scale all entries by a scalar.
    pub fn scale(&mut self, alpha: f64) {
        for v in self.data.iter_mut() {
            *v *= alpha;
        }
    }
}
/// Incomplete LU factorization with zero fill-in (ILU(0)).
///
/// The factorization maintains the same sparsity pattern as the original matrix.
/// Used as a preconditioner for iterative solvers like CG or GMRES.
#[derive(Debug, Clone)]
pub struct Ilu0Preconditioner {
    /// Combined L and U factors stored in CSR format.
    /// L has unit diagonal (not stored), U diagonal is stored.
    /// The storage reuses the sparsity pattern of A.
    pub(super) lu_values: Vec<f64>,
    pub(super) row_ptr: Vec<usize>,
    pub(super) col_indices: Vec<usize>,
    pub(super) n: usize,
}
impl Ilu0Preconditioner {
    /// Compute the ILU(0) factorization of matrix A.
    ///
    /// The factorization uses the same sparsity pattern as A. Elements that
    /// would be fill-in are dropped.
    pub fn new(a: &CsrMatrix) -> Self {
        assert_eq!(a.nrows, a.ncols, "matrix must be square");
        let n = a.nrows;
        let mut lu_values = a.values.clone();
        let row_ptr = a.row_ptr.clone();
        let col_indices = a.col_indices.clone();
        for i in 1..n {
            let row_start = row_ptr[i];
            let row_end = row_ptr[i + 1];
            for p in row_start..row_end {
                let k = col_indices[p];
                if k >= i {
                    break;
                }
                let k_start = row_ptr[k];
                let k_end = row_ptr[k + 1];
                let mut diag_k = 0.0;
                for q in k_start..k_end {
                    if col_indices[q] == k {
                        diag_k = lu_values[q];
                        break;
                    }
                }
                if diag_k.abs() < 1e-60 {
                    continue;
                }
                lu_values[p] /= diag_k;
                let factor = lu_values[p];
                for q in k_start..k_end {
                    let j = col_indices[q];
                    if j <= k {
                        continue;
                    }
                    for s in row_start..row_end {
                        if col_indices[s] == j {
                            lu_values[s] -= factor * lu_values[q];
                            break;
                        }
                    }
                }
            }
        }
        Ilu0Preconditioner {
            lu_values,
            row_ptr,
            col_indices,
            n,
        }
    }
    /// Apply the ILU(0) preconditioner: solve (L U) z = r.
    ///
    /// First forward-substitutes L y = r (L has unit diagonal),
    /// then back-substitutes U z = y.
    pub fn solve(&self, rhs: &[f64]) -> Vec<f64> {
        assert_eq!(rhs.len(), self.n);
        let n = self.n;
        let mut y = rhs.to_vec();
        for i in 0..n {
            let start = self.row_ptr[i];
            let end = self.row_ptr[i + 1];
            for p in start..end {
                let j = self.col_indices[p];
                if j >= i {
                    break;
                }
                y[i] -= self.lu_values[p] * y[j];
            }
        }
        for i in (0..n).rev() {
            let start = self.row_ptr[i];
            let end = self.row_ptr[i + 1];
            let mut diag = 1.0;
            for p in start..end {
                let j = self.col_indices[p];
                if j == i {
                    diag = self.lu_values[p];
                } else if j > i {
                    y[i] -= self.lu_values[p] * y[j];
                }
            }
            if diag.abs() > 1e-60 {
                y[i] /= diag;
            }
        }
        y
    }
}
/// Compressed Sparse Column (CSC) matrix.
///
/// Stores a sparse matrix in CSC format with column pointers, row indices,
/// and values arrays. Efficient for column access and certain direct solvers.
#[derive(Debug, Clone)]
pub struct CscMatrix {
    /// Column pointer array (length = ncols + 1).
    pub col_ptr: Vec<usize>,
    /// Row indices for non-zero entries.
    pub row_indices: Vec<usize>,
    /// Non-zero values.
    pub values: Vec<f64>,
    /// Number of rows.
    pub nrows: usize,
    /// Number of columns.
    pub ncols: usize,
}
impl CscMatrix {
    /// Build a CSC matrix from coordinate (triplet) format.
    ///
    /// Duplicate entries are summed.
    pub fn from_triplets(nrows: usize, ncols: usize, triplets: &[(usize, usize, f64)]) -> Self {
        let mut map: HashMap<(usize, usize), f64> = HashMap::new();
        for &(r, c, v) in triplets {
            assert!(r < nrows, "row index {r} out of bounds for {nrows} rows");
            assert!(c < ncols, "col index {c} out of bounds for {ncols} cols");
            *map.entry((r, c)).or_insert(0.0) += v;
        }
        let mut entries: Vec<((usize, usize), f64)> = map.into_iter().collect();
        entries.sort_by_key(|&((r, c), _)| (c, r));
        let mut col_ptr = vec![0usize; ncols + 1];
        let mut row_indices = Vec::with_capacity(entries.len());
        let mut values = Vec::with_capacity(entries.len());
        for &((r, c), v) in &entries {
            col_ptr[c + 1] += 1;
            row_indices.push(r);
            values.push(v);
        }
        for i in 1..=ncols {
            col_ptr[i] += col_ptr[i - 1];
        }
        Self {
            col_ptr,
            row_indices,
            values,
            nrows,
            ncols,
        }
    }
    /// Get the value at (row, col).
    pub fn get(&self, row: usize, col: usize) -> f64 {
        assert!(row < self.nrows);
        assert!(col < self.ncols);
        let start = self.col_ptr[col];
        let end = self.col_ptr[col + 1];
        for idx in start..end {
            if self.row_indices[idx] == row {
                return self.values[idx];
            }
        }
        0.0
    }
    /// Number of stored non-zero entries.
    pub fn nnz(&self) -> usize {
        self.values.len()
    }
    /// Multiply by a dense vector: y = A * x.
    pub fn mul_vec(&self, x: &[f64]) -> Vec<f64> {
        assert_eq!(x.len(), self.ncols, "vector length must equal ncols");
        let mut y = vec![0.0; self.nrows];
        for (col, &x_col) in x.iter().enumerate().take(self.ncols) {
            let start = self.col_ptr[col];
            let end = self.col_ptr[col + 1];
            for idx in start..end {
                y[self.row_indices[idx]] += self.values[idx] * x_col;
            }
        }
        y
    }
    /// Transpose this CSC matrix, returning a new CSC matrix.
    pub fn transpose(&self) -> CscMatrix {
        let mut triplets = Vec::with_capacity(self.nnz());
        for col in 0..self.ncols {
            let start = self.col_ptr[col];
            let end = self.col_ptr[col + 1];
            for idx in start..end {
                triplets.push((col, self.row_indices[idx], self.values[idx]));
            }
        }
        CscMatrix::from_triplets(self.ncols, self.nrows, &triplets)
    }
    /// Convert to CSR format.
    pub fn to_csr(&self) -> CsrMatrix {
        let mut triplets = Vec::with_capacity(self.nnz());
        for col in 0..self.ncols {
            let start = self.col_ptr[col];
            let end = self.col_ptr[col + 1];
            for idx in start..end {
                triplets.push((self.row_indices[idx], col, self.values[idx]));
            }
        }
        CsrMatrix::from_triplets(self.nrows, self.ncols, &triplets)
    }
}
