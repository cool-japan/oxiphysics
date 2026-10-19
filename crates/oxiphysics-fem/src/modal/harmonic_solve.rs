// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Linear solvers for the real block form of the harmonic-response system.
//!
//! The dynamic stiffness matrix `K − ω² M + iω C` of a harmonic analysis is
//! indefinite above the first natural frequency, so conjugate gradients (which
//! need a positive definite matrix) cannot be used.  This module provides
//!
//! * [`reverse_cuthill_mckee`] — a bandwidth-reducing reordering of the DOF graph,
//! * [`BandLu`] — a banded LU factorisation with partial pivoting (the direct
//!   solver; exact up to rounding for every non-singular system),
//! * [`gmres_ilu0`] — right-preconditioned restarted GMRES with an ILU(0)
//!   preconditioner, the fallback when the band storage of the direct solver
//!   would not fit in memory.
//!
//! All solvers act on a real CSR matrix; the complex harmonic system is embedded
//! as `[[A, −ωC], [ωC, A]]` with interleaved real and imaginary unknowns by
//! [`super::harmonic`].

use crate::error::{Error, Result};
use crate::sparse::{CsrMatrix, Ilu0Preconditioner};

/// Pivots below this value are treated as zero.
///
/// The rows of the matrix are scaled to unit maximum before the factorisation,
/// so the threshold is relative: a matrix whose reciprocal condition number is
/// below about `1e-13` is reported as singular.
const PIVOT_TOL: f64 = 1.0e-13;

/// Largest number of `f64` the band storage of [`BandLu`] may occupy (≈ 200 MB).
pub(super) const BAND_STORAGE_LIMIT: usize = 25_000_000;

// ---------------------------------------------------------------------------
// Reordering
// ---------------------------------------------------------------------------

/// Reverse Cuthill–McKee ordering of an undirected graph.
///
/// `adjacency[i]` lists the neighbours of vertex `i` (self loops are ignored and
/// duplicates tolerated; the graph is symmetrised).  Returns `order` with
/// `order[new] = old`, a permutation of `0..adjacency.len()`.  Every connected
/// component starts from its vertex of minimum degree.
pub(super) fn reverse_cuthill_mckee(adjacency: &[Vec<usize>]) -> Vec<usize> {
    let n = adjacency.len();
    let mut neighbours: Vec<Vec<usize>> = vec![Vec::new(); n];
    for (u, list) in adjacency.iter().enumerate() {
        for &v in list {
            if v != u && v < n {
                neighbours[u].push(v);
                neighbours[v].push(u);
            }
        }
    }
    for list in neighbours.iter_mut() {
        list.sort_unstable();
        list.dedup();
    }
    let degree: Vec<usize> = neighbours.iter().map(Vec::len).collect();

    let mut starts: Vec<usize> = (0..n).collect();
    starts.sort_by_key(|&v| (degree[v], v));

    let mut visited = vec![false; n];
    let mut order: Vec<usize> = Vec::with_capacity(n);
    for &start in &starts {
        if visited[start] {
            continue;
        }
        visited[start] = true;
        let mut head = order.len();
        order.push(start);
        while head < order.len() {
            let u = order[head];
            head += 1;
            let mut next: Vec<usize> = neighbours[u]
                .iter()
                .copied()
                .filter(|&v| !visited[v])
                .collect();
            next.sort_by_key(|&v| (degree[v], v));
            for v in next {
                visited[v] = true;
                order.push(v);
            }
        }
    }
    order.reverse();
    order
}

// ---------------------------------------------------------------------------
// Banded LU
// ---------------------------------------------------------------------------

/// Banded LU factorisation `P A = L U` with row-scaled partial pivoting.
///
/// The matrix is stored in LAPACK general-band layout with `kl` sub-diagonals,
/// `ku` super-diagonals and `kl` extra rows for the fill-in of pivoting (element
/// `(i, j)` lives at `ab[j * ldab + kl + ku + i − j]`).  Cost: `O(n · kl · (kl + ku))`.
pub(super) struct BandLu {
    n: usize,
    kl: usize,
    ku: usize,
    ldab: usize,
    ab: Vec<f64>,
    pivots: Vec<usize>,
    row_scale: Vec<f64>,
}

impl BandLu {
    /// Half bandwidths `(kl, ku)` of the stored pattern of `a`.
    pub(super) fn bandwidths(a: &CsrMatrix) -> (usize, usize) {
        let mut kl = 0usize;
        let mut ku = 0usize;
        for row in 0..a.nrows {
            for &col in &a.col_indices[a.row_ptr[row]..a.row_ptr[row + 1]] {
                if row > col {
                    kl = kl.max(row - col);
                } else {
                    ku = ku.max(col - row);
                }
            }
        }
        (kl, ku)
    }

    /// Number of `f64` needed by the band storage of an `n x n` matrix.
    pub(super) fn storage_len(n: usize, kl: usize, ku: usize) -> usize {
        n.saturating_mul(2 * kl + ku + 1)
    }

    /// Factor the square matrix `a`.
    ///
    /// # Errors
    ///
    /// * [`Error::InvalidInput`] if `a` is not square or holds a non-finite entry.
    /// * [`Error::SingularSystem`] if a row is zero or a pivot falls below
    ///   [`PIVOT_TOL`] (relative to the row-scaled matrix).
    pub(super) fn factor(a: &CsrMatrix) -> Result<Self> {
        if a.nrows != a.ncols {
            return Err(Error::InvalidInput(format!(
                "banded LU needs a square matrix, got {} x {}",
                a.nrows, a.ncols
            )));
        }
        let n = a.nrows;
        let (kl, ku) = Self::bandwidths(a);
        let kv = kl + ku;
        let ldab = 2 * kl + ku + 1;

        let mut row_scale = vec![1.0_f64; n];
        for (row, scale) in row_scale.iter_mut().enumerate() {
            let range = a.row_ptr[row]..a.row_ptr[row + 1];
            let mut largest = 0.0_f64;
            for &value in &a.values[range] {
                if !value.is_finite() {
                    return Err(Error::InvalidInput(
                        "the system matrix holds a non-finite entry".to_string(),
                    ));
                }
                largest = largest.max(value.abs());
            }
            if largest == 0.0 {
                return Err(Error::SingularSystem);
            }
            *scale = 1.0 / largest;
        }

        let mut ab = vec![0.0_f64; n * ldab];
        for row in 0..n {
            for idx in a.row_ptr[row]..a.row_ptr[row + 1] {
                let col = a.col_indices[idx];
                ab[col * ldab + (kv + row - col)] += a.values[idx] * row_scale[row];
            }
        }

        let mut pivots: Vec<usize> = (0..n).collect();
        let mut last_col = 0usize;
        for (j, pivot_slot) in pivots.iter_mut().enumerate() {
            let km = kl.min(n - 1 - j);
            let diag = j * ldab + kv;
            let mut jp = 0usize;
            let mut best = ab[diag].abs();
            for t in 1..=km {
                let candidate = ab[diag + t].abs();
                if candidate > best {
                    best = candidate;
                    jp = t;
                }
            }
            if best.is_nan() || best <= PIVOT_TOL {
                return Err(Error::SingularSystem);
            }
            *pivot_slot = j + jp;
            last_col = last_col.max((j + ku + jp).min(n - 1));
            if jp != 0 {
                for col in j..=last_col {
                    let upper = col * ldab + (kv + j - col);
                    ab.swap(upper, upper + jp);
                }
            }
            if km > 0 {
                let pivot = ab[diag];
                for t in 1..=km {
                    ab[diag + t] /= pivot;
                }
                for col in (j + 1)..=last_col {
                    let upper = col * ldab + (kv + j - col);
                    let u = ab[upper];
                    if u != 0.0 {
                        for t in 1..=km {
                            ab[upper + t] -= ab[diag + t] * u;
                        }
                    }
                }
            }
        }
        Ok(Self {
            n,
            kl,
            ku,
            ldab,
            ab,
            pivots,
            row_scale,
        })
    }

    /// Solve `A x = rhs` with the stored factors.
    pub(super) fn solve(&self, rhs: &[f64]) -> Vec<f64> {
        let n = self.n;
        let kv = self.kl + self.ku;
        let mut b: Vec<f64> = rhs
            .iter()
            .zip(self.row_scale.iter())
            .map(|(value, scale)| value * scale)
            .collect();
        // Forward substitution: apply the row interchanges and the unit lower factor.
        for j in 0..n.saturating_sub(1) {
            let lm = self.kl.min(n - 1 - j);
            let pivot_row = self.pivots[j];
            if pivot_row != j {
                b.swap(pivot_row, j);
            }
            let bj = b[j];
            if bj != 0.0 {
                let column = j * self.ldab + kv;
                for t in 1..=lm {
                    b[j + t] -= self.ab[column + t] * bj;
                }
            }
        }
        // Back substitution with the banded upper factor (bandwidth kl + ku).
        for j in (0..n).rev() {
            let column = j * self.ldab + kv;
            b[j] /= self.ab[column];
            let bj = b[j];
            if bj != 0.0 {
                let first = j.saturating_sub(kv);
                for (i, bi) in b.iter_mut().enumerate().take(j).skip(first) {
                    *bi -= self.ab[column - (j - i)] * bj;
                }
            }
        }
        b
    }
}

// ---------------------------------------------------------------------------
// ILU(0)-preconditioned GMRES
// ---------------------------------------------------------------------------

/// Outcome of an iterative solve: the iterate and its *true* residual.
pub(super) struct KrylovOutcome {
    /// Last iterate.
    pub x: Vec<f64>,
    /// Number of Arnoldi steps performed.
    pub iterations: usize,
    /// True relative residual `‖b − A x‖ / ‖b‖` of [`Self::x`].
    pub relative_residual: f64,
}

fn norm2(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

fn axpy(alpha: f64, x: &[f64], y: &mut [f64]) {
    for (yi, xi) in y.iter_mut().zip(x.iter()) {
        *yi += alpha * xi;
    }
}

/// Copy of `a` whose diagonal is present and bounded away from zero, used only
/// to build the ILU(0) preconditioner (a zero or tiny pivot would otherwise make
/// the incomplete factorisation meaningless; the iteration still solves the
/// original system).
fn with_safeguarded_diagonal(a: &CsrMatrix) -> CsrMatrix {
    const RELATIVE_FLOOR: f64 = 1.0e-2;
    let mut triplets: Vec<(usize, usize, f64)> = Vec::with_capacity(a.nnz() + a.nrows);
    for row in 0..a.nrows {
        let range = a.row_ptr[row]..a.row_ptr[row + 1];
        let mut largest = 0.0_f64;
        let mut diagonal = 0.0_f64;
        for idx in range {
            let col = a.col_indices[idx];
            let value = a.values[idx];
            largest = largest.max(value.abs());
            if col == row {
                diagonal = value;
            }
            triplets.push((row, col, value));
        }
        let floor = RELATIVE_FLOOR * largest;
        if floor > 0.0 && diagonal.abs() < floor {
            let target = if diagonal < 0.0 { -floor } else { floor };
            triplets.push((row, row, target - diagonal));
        }
    }
    CsrMatrix::from_triplets(a.nrows, a.ncols, &triplets)
}

/// Right-preconditioned restarted GMRES(`restart`) with an ILU(0) preconditioner.
///
/// Solves `A x = b` for a general (non-symmetric, possibly indefinite) matrix.
/// A cycle is stopped when the Givens estimate of the preconditioned residual
/// falls below `tol · ‖b‖`; the residual reported in the outcome is always the
/// **true** residual `‖b − A x‖ / ‖b‖` recomputed from the iterate, so a
/// stagnating or broken-down iteration cannot masquerade as converged.
pub(super) fn gmres_ilu0(
    a: &CsrMatrix,
    b: &[f64],
    restart: usize,
    max_cycles: usize,
    tol: f64,
) -> KrylovOutcome {
    let n = b.len();
    let b_norm = norm2(b);
    let mut x = vec![0.0_f64; n];
    if n == 0 || b_norm == 0.0 {
        return KrylovOutcome {
            x,
            iterations: 0,
            relative_residual: 0.0,
        };
    }
    let precond = Ilu0Preconditioner::new(&with_safeguarded_diagonal(a));
    let m = restart.clamp(1, n);
    let mut iterations = 0usize;

    let true_residual = |x: &[f64]| -> (Vec<f64>, f64) {
        let ax = a.mul_vec(x);
        let r: Vec<f64> = b.iter().zip(ax.iter()).map(|(bi, ai)| bi - ai).collect();
        let rel = norm2(&r) / b_norm;
        (r, rel)
    };

    let (mut r, mut rel) = true_residual(&x);
    for _cycle in 0..max_cycles {
        if rel.is_nan() || rel <= tol {
            break;
        }
        let beta = rel * b_norm;
        let mut basis: Vec<Vec<f64>> = Vec::with_capacity(m + 1);
        basis.push(r.iter().map(|v| v / beta).collect());
        let mut h = vec![vec![0.0_f64; m]; m + 1];
        let mut cs = vec![0.0_f64; m];
        let mut sn = vec![0.0_f64; m];
        let mut g = vec![0.0_f64; m + 1];
        g[0] = beta;
        let mut used = 0usize;
        for j in 0..m {
            let z = precond.solve(&basis[j]);
            let mut w = a.mul_vec(&z);
            for i in 0..=j {
                let hij = dot(&w, &basis[i]);
                h[i][j] = hij;
                axpy(-hij, &basis[i], &mut w);
            }
            let w_norm = norm2(&w);
            h[j + 1][j] = w_norm;
            for i in 0..j {
                let upper = cs[i] * h[i][j] + sn[i] * h[i + 1][j];
                h[i + 1][j] = -sn[i] * h[i][j] + cs[i] * h[i + 1][j];
                h[i][j] = upper;
            }
            let radius = h[j][j].hypot(h[j + 1][j]);
            if radius > 0.0 {
                cs[j] = h[j][j] / radius;
                sn[j] = h[j + 1][j] / radius;
            } else {
                cs[j] = 1.0;
                sn[j] = 0.0;
            }
            h[j][j] = radius;
            h[j + 1][j] = 0.0;
            g[j + 1] = -sn[j] * g[j];
            g[j] *= cs[j];
            iterations += 1;
            used = j + 1;
            if !w_norm.is_finite() || w_norm <= 1.0e-14 * beta {
                break;
            }
            if g[j + 1].abs() <= tol * b_norm {
                break;
            }
            basis.push(w.iter().map(|v| v / w_norm).collect());
        }
        let mut y = vec![0.0_f64; used];
        for i in (0..used).rev() {
            let mut sum = g[i];
            for k in (i + 1)..used {
                sum -= h[i][k] * y[k];
            }
            y[i] = if h[i][i] != 0.0 { sum / h[i][i] } else { 0.0 };
        }
        let mut combination = vec![0.0_f64; n];
        for (yi, vi) in y.iter().zip(basis.iter()) {
            axpy(*yi, vi, &mut combination);
        }
        let correction = precond.solve(&combination);
        axpy(1.0, &correction, &mut x);
        let (new_r, new_rel) = true_residual(&x);
        r = new_r;
        rel = new_rel;
    }
    KrylovOutcome {
        x,
        iterations,
        relative_residual: rel,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Deterministic pseudo-random numbers in `[-1, 1)` (xorshift64*).
    struct Lcg(u64);

    impl Lcg {
        fn next_f64(&mut self) -> f64 {
            self.0 ^= self.0 >> 12;
            self.0 ^= self.0 << 25;
            self.0 ^= self.0 >> 27;
            let bits = self.0.wrapping_mul(0x2545_F491_4F6C_DD1D) >> 11;
            2.0 * (bits as f64 / (1u64 << 53) as f64) - 1.0
        }
    }

    /// Dense Gaussian elimination with full row pivoting (reference solver).
    fn dense_solve(a: &[Vec<f64>], b: &[f64]) -> Option<Vec<f64>> {
        let n = b.len();
        let mut m: Vec<Vec<f64>> = a.to_vec();
        let mut rhs = b.to_vec();
        for col in 0..n {
            let mut p = col;
            for row in (col + 1)..n {
                if m[row][col].abs() > m[p][col].abs() {
                    p = row;
                }
            }
            if m[p][col].abs() < 1e-300 {
                return None;
            }
            m.swap(col, p);
            rhs.swap(col, p);
            for row in (col + 1)..n {
                let f = m[row][col] / m[col][col];
                let pivot_values = m[col].clone();
                for (entry, upper) in m[row].iter_mut().zip(pivot_values.iter()).skip(col) {
                    *entry -= f * upper;
                }
                let v = rhs[col];
                rhs[row] -= f * v;
            }
        }
        let mut x = vec![0.0; n];
        for row in (0..n).rev() {
            let mut s = rhs[row];
            for k in (row + 1)..n {
                s -= m[row][k] * x[k];
            }
            x[row] = s / m[row][row];
        }
        Some(x)
    }

    fn to_csr(dense: &[Vec<f64>]) -> CsrMatrix {
        let n = dense.len();
        let mut triplets = Vec::new();
        for (i, row) in dense.iter().enumerate() {
            for (j, &v) in row.iter().enumerate() {
                if v != 0.0 {
                    triplets.push((i, j, v));
                }
            }
        }
        CsrMatrix::from_triplets(n, n, &triplets)
    }

    /// Random non-symmetric banded matrix with a *zero* diagonal in several rows,
    /// so the factorisation cannot succeed without pivoting.
    fn random_banded(n: usize, kl: usize, ku: usize, rng: &mut Lcg) -> Vec<Vec<f64>> {
        let mut a = vec![vec![0.0; n]; n];
        for (i, row) in a.iter_mut().enumerate() {
            let lo = i.saturating_sub(kl);
            let hi = (i + ku).min(n - 1);
            for entry in row.iter_mut().take(hi + 1).skip(lo) {
                *entry = rng.next_f64();
            }
            if i % 3 == 0 {
                row[i] = 0.0;
            }
        }
        a
    }

    #[test]
    fn band_lu_matches_dense_reference_on_pivoting_matrices() {
        let mut rng = Lcg(0x9E37_79B9_7F4A_7C15);
        for &(n, kl, ku) in &[
            (1usize, 0usize, 0usize),
            (2, 1, 1),
            (7, 2, 1),
            (12, 1, 3),
            (30, 4, 4),
            (25, 0, 3),
            (25, 3, 0),
        ] {
            let dense = random_banded(n, kl, ku, &mut rng);
            let rhs: Vec<f64> = (0..n).map(|_| rng.next_f64()).collect();
            let reference = match dense_solve(&dense, &rhs) {
                Some(x) => x,
                None => continue,
            };
            let a = to_csr(&dense);
            let lu = BandLu::factor(&a).expect("random matrix is non-singular");
            let x = lu.solve(&rhs);
            for (got, want) in x.iter().zip(reference.iter()) {
                assert!(
                    (got - want).abs() <= 1e-9 * (1.0 + want.abs()),
                    "n={n} kl={kl} ku={ku}: {got} vs {want}"
                );
            }
            // The residual of the banded solve must be at rounding level.
            let ax = a.mul_vec(&x);
            let resid = ax
                .iter()
                .zip(rhs.iter())
                .map(|(l, r)| (l - r).abs())
                .fold(0.0_f64, f64::max);
            assert!(resid < 1e-9, "n={n} residual {resid}");
        }
    }

    #[test]
    fn band_lu_detects_singular_matrices() {
        // Two identical rows.
        let dense = vec![
            vec![1.0, 2.0, 0.0],
            vec![1.0, 2.0, 0.0],
            vec![0.0, 1.0, 3.0],
        ];
        assert!(matches!(
            BandLu::factor(&to_csr(&dense)),
            Err(Error::SingularSystem)
        ));
        // A zero row.
        let zero_row = vec![vec![1.0, 0.0], vec![0.0, 0.0]];
        assert!(matches!(
            BandLu::factor(&to_csr(&zero_row)),
            Err(Error::SingularSystem)
        ));
        // Numerically singular: second row is a multiple of the first up to 1e-16.
        let nearly = vec![vec![1.0, 1.0], vec![1.0, 1.0 + 1.0e-16]];
        assert!(matches!(
            BandLu::factor(&to_csr(&nearly)),
            Err(Error::SingularSystem)
        ));
    }

    #[test]
    fn band_lu_handles_badly_scaled_rows() {
        // Rows differ by 20 orders of magnitude; row scaling keeps the pivot test meaningful.
        let dense = vec![
            vec![2.0e10, -1.0e10, 0.0],
            vec![-1.0e-10, 2.0e-10, -1.0e-10],
            vec![0.0, -1.0, 2.0],
        ];
        let rhs = vec![1.0e10, 0.0, 1.0];
        let reference = dense_solve(&dense, &rhs).expect("non-singular");
        let lu = BandLu::factor(&to_csr(&dense)).expect("non-singular after scaling");
        let x = lu.solve(&rhs);
        for (got, want) in x.iter().zip(reference.iter()) {
            assert!(
                (got - want).abs() <= 1e-9 * (1.0 + want.abs()),
                "{got} vs {want}"
            );
        }
    }

    #[test]
    fn band_lu_rejects_non_square_and_non_finite_input() {
        let rect = CsrMatrix::from_triplets(2, 3, &[(0, 0, 1.0), (1, 1, 1.0)]);
        assert!(matches!(BandLu::factor(&rect), Err(Error::InvalidInput(_))));
        let bad = CsrMatrix::from_triplets(2, 2, &[(0, 0, f64::NAN), (1, 1, 1.0)]);
        assert!(matches!(BandLu::factor(&bad), Err(Error::InvalidInput(_))));
    }

    #[test]
    fn bandwidths_and_storage_of_a_tridiagonal_pattern() {
        let a = to_csr(&[
            vec![2.0, -1.0, 0.0, 0.0],
            vec![-1.0, 2.0, -1.0, 0.0],
            vec![0.0, -1.0, 2.0, -1.0],
            vec![0.0, 0.0, -1.0, 2.0],
        ]);
        assert_eq!(BandLu::bandwidths(&a), (1, 1));
        assert_eq!(BandLu::storage_len(4, 1, 1), 4 * 4);
    }

    /// Bandwidth `max |new(i) − new(j)|` over the edges of the graph.
    fn bandwidth_of(adjacency: &[Vec<usize>], order: &[usize]) -> usize {
        let mut position = vec![0usize; order.len()];
        for (new, &old) in order.iter().enumerate() {
            position[old] = new;
        }
        let mut width = 0usize;
        for (u, list) in adjacency.iter().enumerate() {
            for &v in list {
                width = width.max(position[u].abs_diff(position[v]));
            }
        }
        width
    }

    #[test]
    fn rcm_restores_unit_bandwidth_of_a_scrambled_chain() {
        // A path graph whose vertices are numbered in a scrambled order.
        let n = 41usize;
        let label: Vec<usize> = (0..n).map(|i| (i * 17) % n).collect();
        let mut adjacency = vec![Vec::new(); n];
        for i in 0..n - 1 {
            adjacency[label[i]].push(label[i + 1]);
            adjacency[label[i + 1]].push(label[i]);
        }
        let identity: Vec<usize> = (0..n).collect();
        assert!(bandwidth_of(&adjacency, &identity) > 10);
        let order = reverse_cuthill_mckee(&adjacency);
        let mut sorted = order.clone();
        sorted.sort_unstable();
        assert_eq!(sorted, identity, "RCM must return a permutation");
        assert_eq!(bandwidth_of(&adjacency, &order), 1);
    }

    #[test]
    fn rcm_orders_every_vertex_of_a_disconnected_graph() {
        // Two components plus an isolated vertex; duplicate and self edges are tolerated.
        let adjacency = vec![vec![1, 1, 0], vec![0, 2], vec![1], vec![4], vec![3], vec![]];
        let order = reverse_cuthill_mckee(&adjacency);
        let mut sorted = order.clone();
        sorted.sort_unstable();
        assert_eq!(sorted, vec![0, 1, 2, 3, 4, 5]);
        assert!(bandwidth_of(&adjacency, &order) <= 1);
        assert!(reverse_cuthill_mckee(&[]).is_empty());
    }

    #[test]
    fn rcm_reduces_the_bandwidth_of_a_scrambled_grid() {
        let (nx, ny) = (9usize, 6usize);
        let n = nx * ny;
        let label: Vec<usize> = (0..n).map(|i| (i * 29) % n).collect();
        let mut adjacency = vec![Vec::new(); n];
        for j in 0..ny {
            for i in 0..nx {
                let u = label[j * nx + i];
                if i + 1 < nx {
                    adjacency[u].push(label[j * nx + i + 1]);
                }
                if j + 1 < ny {
                    adjacency[u].push(label[(j + 1) * nx + i]);
                }
            }
        }
        let identity: Vec<usize> = (0..n).collect();
        let scrambled = bandwidth_of(&adjacency, &identity);
        let order = reverse_cuthill_mckee(&adjacency);
        let reordered = bandwidth_of(&adjacency, &order);
        assert!(
            reordered <= 2 * nx.min(ny) + 2 && reordered * 2 < scrambled,
            "bandwidth {scrambled} -> {reordered}"
        );
    }

    /// Damped dynamic-stiffness-like chain `[[A, -B], [B, A]]` with interleaved unknowns.
    fn interleaved_chain(n: usize, a_diag: f64, coupling: f64, damping: f64) -> CsrMatrix {
        let mut triplets = Vec::new();
        for i in 0..n {
            let mut push_block = |row: usize, col: usize, a: f64, b: f64| {
                triplets.push((2 * row, 2 * col, a));
                triplets.push((2 * row + 1, 2 * col + 1, a));
                triplets.push((2 * row, 2 * col + 1, -b));
                triplets.push((2 * row + 1, 2 * col, b));
            };
            push_block(i, i, a_diag, damping);
            if i + 1 < n {
                push_block(i, i + 1, coupling, 0.1 * damping);
                push_block(i + 1, i, coupling, 0.1 * damping);
            }
        }
        CsrMatrix::from_triplets(2 * n, 2 * n, &triplets)
    }

    #[test]
    fn gmres_ilu0_matches_the_direct_solution_of_a_damped_chain() {
        let n = 60usize;
        // Diagonally dominant, so ILU(0)-GMRES converges rapidly.
        let a = interleaved_chain(n, 4.0, -1.0, 0.3);
        let rhs: Vec<f64> = (0..2 * n).map(|i| ((i % 7) as f64) - 3.0).collect();
        let direct = BandLu::factor(&a).expect("non-singular").solve(&rhs);
        let outcome = gmres_ilu0(&a, &rhs, 30, 50, 1.0e-12);
        assert!(
            outcome.relative_residual <= 1.0e-10,
            "residual {}",
            outcome.relative_residual
        );
        assert!(outcome.iterations > 0);
        for (got, want) in outcome.x.iter().zip(direct.iter()) {
            assert!(
                (got - want).abs() <= 1e-8 * (1.0 + want.abs()),
                "{got} vs {want}"
            );
        }
    }

    #[test]
    fn gmres_ilu0_solves_an_indefinite_system_and_reports_true_residual() {
        // Indefinite: the diagonal alternates in sign (zero-free), moderate size.
        let n = 24usize;
        let mut triplets = Vec::new();
        for i in 0..n {
            let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
            triplets.push((i, i, 3.0 * sign));
            if i + 1 < n {
                triplets.push((i, i + 1, -1.0));
                triplets.push((i + 1, i, -0.5));
            }
        }
        let a = CsrMatrix::from_triplets(n, n, &triplets);
        let rhs: Vec<f64> = (0..n).map(|i| 1.0 + 0.1 * i as f64).collect();
        let outcome = gmres_ilu0(&a, &rhs, n, 20, 1.0e-12);
        let ax = a.mul_vec(&outcome.x);
        let true_rel = norm2(
            &ax.iter()
                .zip(rhs.iter())
                .map(|(l, r)| l - r)
                .collect::<Vec<f64>>(),
        ) / norm2(&rhs);
        assert!(true_rel < 1.0e-10, "true relative residual {true_rel}");
        assert!(
            (outcome.relative_residual - true_rel).abs() <= 1e-12 + 1e-6 * true_rel,
            "reported {} vs recomputed {true_rel}",
            outcome.relative_residual
        );
    }

    #[test]
    fn gmres_ilu0_reports_failure_instead_of_claiming_convergence() {
        // One cycle of two inner steps cannot solve a 60-unknown indefinite system.
        let n = 60usize;
        let mut triplets = Vec::new();
        for i in 0..n {
            let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
            triplets.push((i, i, 0.2 * sign));
            if i + 1 < n {
                triplets.push((i, i + 1, 1.0));
                triplets.push((i + 1, i, -1.0));
            }
            if i + 5 < n {
                triplets.push((i, i + 5, 0.7));
            }
        }
        let a = CsrMatrix::from_triplets(n, n, &triplets);
        let rhs: Vec<f64> = (0..n).map(|i| (i as f64).sin() + 2.0).collect();
        let outcome = gmres_ilu0(&a, &rhs, 2, 1, 1.0e-12);
        assert!(
            outcome.relative_residual > 1.0e-6,
            "two Arnoldi steps cannot converge, got {}",
            outcome.relative_residual
        );
        assert_eq!(outcome.iterations, 2);
    }

    #[test]
    fn gmres_ilu0_zero_rhs_and_empty_system() {
        let a = CsrMatrix::from_triplets(2, 2, &[(0, 0, 1.0), (1, 1, 1.0)]);
        let out = gmres_ilu0(&a, &[0.0, 0.0], 4, 4, 1e-12);
        assert_eq!(out.x, vec![0.0, 0.0]);
        assert_eq!(out.relative_residual, 0.0);
        let empty = CsrMatrix::new(0, 0);
        let out = gmres_ilu0(&empty, &[], 4, 4, 1e-12);
        assert!(out.x.is_empty());
    }

    #[test]
    fn safeguarded_diagonal_inserts_missing_and_tiny_pivots_only() {
        let a = CsrMatrix::from_triplets(
            3,
            3,
            &[
                (0, 1, 2.0),
                (1, 1, 1.0e-9),
                (1, 2, 3.0),
                (2, 2, 5.0),
                (2, 0, 1.0),
            ],
        );
        let s = with_safeguarded_diagonal(&a);
        assert!(
            (s.get(0, 0) - 2.0e-2).abs() < 1e-15,
            "missing diagonal inserted"
        );
        assert!((s.get(1, 1) - 3.0e-2).abs() < 1e-15, "tiny diagonal lifted");
        assert!(
            (s.get(2, 2) - 5.0).abs() < 1e-15,
            "healthy diagonal untouched"
        );
        assert!((s.get(0, 1) - 2.0).abs() < 1e-15);
    }
}
