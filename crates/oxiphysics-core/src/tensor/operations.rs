//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use super::types::*;

/// Einstein summation helper for common rank-2 tensor operations.
///
/// Supported notation strings:
/// - `"ij,jk->ik"` : matrix–matrix product
/// - `"ij,ij->"` : Frobenius inner product (returns scalar in a 1-element vec)
/// - `"ij->ji"` : transpose
/// - `"ii->"` : trace (returns scalar in a 1-element vec)
#[allow(dead_code)]
pub fn einsum_2d(notation: &str, a: &[Vec<f64>], b: Option<&[Vec<f64>]>) -> Vec<Vec<f64>> {
    let notation = notation.trim();
    match notation {
        "ij,jk->ik" => {
            let b = b.expect("einsum 'ij,jk->ik' requires two operands");
            let m = a.len();
            let k = a[0].len();
            let n = b[0].len();
            let mut c = vec![vec![0.0; n]; m];
            for i in 0..m {
                for kk in 0..k {
                    for j in 0..n {
                        c[i][j] += a[i][kk] * b[kk][j];
                    }
                }
            }
            c
        }
        "ij,ij->" => {
            let b = b.expect("einsum 'ij,ij->' requires two operands");
            let s: f64 = a
                .iter()
                .zip(b.iter())
                .flat_map(|(ar, br)| ar.iter().zip(br.iter()).map(|(&x, &y)| x * y))
                .sum();
            vec![vec![s]]
        }
        "ij->ji" => {
            let m = a.len();
            let n = a[0].len();
            let mut out = vec![vec![0.0; m]; n];
            for i in 0..m {
                for j in 0..n {
                    out[j][i] = a[i][j];
                }
            }
            out
        }
        "ii->" => {
            let s: f64 = a.iter().enumerate().map(|(i, row)| row[i]).sum();
            vec![vec![s]]
        }
        _ => panic!("einsum_2d: unsupported notation '{notation}'"),
    }
}
/// Helper: multiply m×k matrix by k×n matrix.
#[allow(dead_code)]
pub(super) fn matmul(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let m = a.len();
    let k = a[0].len();
    let n = b[0].len();
    let mut c = vec![vec![0.0; n]; m];
    for i in 0..m {
        for p in 0..k {
            for j in 0..n {
                c[i][j] += a[i][p] * b[p][j];
            }
        }
    }
    c
}
/// Helper: element-wise (Hadamard) product of two matrices of the same shape.
#[allow(dead_code)]
fn hadamard(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    a.iter()
        .zip(b.iter())
        .map(|(ra, rb)| ra.iter().zip(rb.iter()).map(|(&x, &y)| x * y).collect())
        .collect()
}
/// Helper: transpose a matrix.
#[allow(dead_code)]
pub(super) fn transpose(a: &[Vec<f64>]) -> Vec<Vec<f64>> {
    if a.is_empty() {
        return vec![];
    }
    let m = a.len();
    let n = a[0].len();
    let mut t = vec![vec![0.0; m]; n];
    for i in 0..m {
        for j in 0..n {
            t[j][i] = a[i][j];
        }
    }
    t
}
/// Helper: solve a small linear system A x = b via Gaussian elimination (dense).
#[allow(dead_code)]
fn solve_ls(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let r = a.len();
    let m = b[0].len();
    let mut aug: Vec<Vec<f64>> = a
        .iter()
        .enumerate()
        .map(|(i, row)| {
            let mut r = row.clone();
            r.extend_from_slice(&b[i]);
            r
        })
        .collect();
    for col in 0..r {
        let mut pivot = col;
        for row in (col + 1)..r {
            if aug[row][col].abs() > aug[pivot][col].abs() {
                pivot = row;
            }
        }
        aug.swap(col, pivot);
        let d = aug[col][col];
        if d.abs() < 1e-14 {
            continue;
        }
        for j in col..r + m {
            aug[col][j] /= d;
        }
        for row in 0..r {
            if row == col {
                continue;
            }
            let factor = aug[row][col];
            for j in col..r + m {
                aug[row][j] -= factor * aug[col][j];
            }
        }
    }
    let mut x = vec![vec![0.0; m]; r];
    for i in 0..r {
        for j in 0..m {
            x[i][j] = aug[i][r + j];
        }
    }
    x
}
/// Alternating Least Squares (ALS) for CP decomposition of a 3-way tensor.
///
/// Decomposes `tensor` (shape n0 × n1 × n2) into `rank` components.
/// Returns the factor matrices A (n0×rank), B (n1×rank), C (n2×rank) and
/// normalisation weights λ.
///
/// # Arguments
/// * `tensor`   – reference to a `DenseTensor` with rank 3
/// * `rank`     – number of CP components
/// * `max_iter` – maximum ALS iterations
/// * `tol`      – relative reconstruction error tolerance for convergence
#[allow(dead_code)]
pub fn cp_als(tensor: &DenseTensor, rank: usize, max_iter: usize, tol: f64) -> CpDecomposition {
    assert_eq!(tensor.shape.len(), 3, "cp_als requires a rank-3 tensor");
    let n0 = tensor.shape[0];
    let n1 = tensor.shape[1];
    let n2 = tensor.shape[2];
    let init_val = 1.0 / rank as f64;
    let mut a: Vec<Vec<f64>> = (0..n0)
        .map(|i| {
            (0..rank)
                .map(|r| ((i + r) as f64 + 1.0) * init_val)
                .collect()
        })
        .collect();
    let mut b: Vec<Vec<f64>> = (0..n1)
        .map(|j| (0..rank).map(|r| ((j + r + 1) as f64) * init_val).collect())
        .collect();
    let mut c: Vec<Vec<f64>> = (0..n2)
        .map(|k| (0..rank).map(|r| ((k + r + 2) as f64) * init_val).collect())
        .collect();
    let x0 = tensor.mode_n_unfold(0);
    let x1 = tensor.mode_n_unfold(1);
    let x2 = tensor.mode_n_unfold(2);
    let mut prev_err = f64::MAX;
    for _iter in 0..max_iter {
        let c_tc = gram(&c);
        let b_tb = gram(&b);
        let v_cb = hadamard(&c_tc, &b_tb);
        let khatri_rao_cb = khatri_rao(&c, &b);
        let rhs_a = matmul(&x0, &khatri_rao_cb);
        let rhs_at = transpose(&rhs_a);
        let sol_a = solve_ls(&v_cb, &rhs_at);
        a = transpose(&sol_a);
        let a_ta = gram(&a);
        let v_ca = hadamard(&c_tc, &a_ta);
        let khatri_rao_ca = khatri_rao(&c, &a);
        let rhs_b = matmul(&x1, &khatri_rao_ca);
        let rhs_bt = transpose(&rhs_b);
        let sol_b = solve_ls(&v_ca, &rhs_bt);
        b = transpose(&sol_b);
        let b_tb2 = gram(&b);
        let a_ta2 = gram(&a);
        let v_ba = hadamard(&b_tb2, &a_ta2);
        let khatri_rao_ba = khatri_rao(&b, &a);
        let rhs_c = matmul(&x2, &khatri_rao_ba);
        let rhs_ct = transpose(&rhs_c);
        let sol_c = solve_ls(&v_ba, &rhs_ct);
        c = transpose(&sol_c);
        let mut lambdas = vec![1.0f64; rank];
        for r in 0..rank {
            let na: f64 = a.iter().map(|row| row[r] * row[r]).sum::<f64>().sqrt();
            let nb: f64 = b.iter().map(|row| row[r] * row[r]).sum::<f64>().sqrt();
            let nc: f64 = c.iter().map(|row| row[r] * row[r]).sum::<f64>().sqrt();
            let lam = na * nb * nc;
            lambdas[r] = lam;
            if na > 1e-14 {
                for row in &mut a {
                    row[r] /= na;
                }
            }
            if nb > 1e-14 {
                for row in &mut b {
                    row[r] /= nb;
                }
            }
            if nc > 1e-14 {
                for row in &mut c {
                    row[r] /= nc;
                }
            }
        }
        let err = cp_reconstruction_error(tensor, &a, &b, &c, &lambdas);
        if (prev_err - err).abs() / (prev_err.abs() + 1e-14) < tol {
            return CpDecomposition { a, b, c, lambdas };
        }
        prev_err = err;
    }
    let lambdas = vec![1.0f64; rank];
    CpDecomposition { a, b, c, lambdas }
}
/// Compute Gram matrix A^T A of a factor matrix A (n×r) → (r×r).
#[allow(dead_code)]
fn gram(a: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let r = a[0].len();
    let mut g = vec![vec![0.0; r]; r];
    for row in a {
        for i in 0..r {
            for j in 0..r {
                g[i][j] += row[i] * row[j];
            }
        }
    }
    g
}
/// Khatri-Rao product (column-wise Kronecker product) of A (m×r) and B (n×r)
/// → (m*n)×r.
#[allow(dead_code)]
fn khatri_rao(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let m = a.len();
    let n = b.len();
    let r = a[0].len();
    let mut out = vec![vec![0.0; r]; m * n];
    for i in 0..m {
        for j in 0..n {
            for rr in 0..r {
                out[i * n + j][rr] = a[i][rr] * b[j][rr];
            }
        }
    }
    out
}
/// Frobenius norm of the reconstruction error for a CP decomposition.
#[allow(dead_code)]
fn cp_reconstruction_error(
    tensor: &DenseTensor,
    a: &[Vec<f64>],
    b: &[Vec<f64>],
    c: &[Vec<f64>],
    lambdas: &[f64],
) -> f64 {
    let n0 = tensor.shape[0];
    let n1 = tensor.shape[1];
    let n2 = tensor.shape[2];
    let rank = lambdas.len();
    let mut err = 0.0;
    for i in 0..n0 {
        for j in 0..n1 {
            for k in 0..n2 {
                let approx: f64 = (0..rank)
                    .map(|r| lambdas[r] * a[i][r] * b[j][r] * c[k][r])
                    .sum();
                let diff = tensor.get(&[i, j, k]) - approx;
                err += diff * diff;
            }
        }
    }
    err.sqrt()
}
/// Helper: extract leading `k` left singular vectors of a matrix via power
/// iteration on A A^T with Gram-Schmidt orthogonalisation.  Returns an n×k matrix.
#[allow(dead_code)]
pub(super) fn truncated_svd_left(a: &[Vec<f64>], k: usize) -> Vec<Vec<f64>> {
    let n = a.len();
    let m = a[0].len();
    let mut aat = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            for p in 0..m {
                aat[i][j] += a[i][p] * a[j][p];
            }
        }
    }
    let k = k.min(n);
    let mut result = vec![vec![0.0f64; k]; n];
    let mut found: Vec<Vec<f64>> = Vec::new();
    for r in 0..k {
        let mut v: Vec<f64> = (0..n).map(|i| if i == r { 1.0 } else { 0.0 }).collect();
        for fv in &found {
            let dot: f64 = v.iter().zip(fv.iter()).map(|(a, b)| a * b).sum();
            for i in 0..n {
                v[i] -= dot * fv[i];
            }
        }
        let norm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        if norm < 1e-14 {
            v = (0..n)
                .map(|i| if i == (r + 1) % n { 1.0 } else { 0.0 })
                .collect();
        } else {
            for vi in &mut v {
                *vi /= norm;
            }
        }
        for _ in 0..300 {
            let mut w = vec![0.0; n];
            for i in 0..n {
                for j in 0..n {
                    w[i] += aat[i][j] * v[j];
                }
            }
            for fv in &found {
                let dot: f64 = w.iter().zip(fv.iter()).map(|(a, b)| a * b).sum();
                for i in 0..n {
                    w[i] -= dot * fv[i];
                }
            }
            let norm = w.iter().map(|x| x * x).sum::<f64>().sqrt();
            if norm < 1e-14 {
                break;
            }
            let v_new: Vec<f64> = w.iter().map(|x| x / norm).collect();
            let diff: f64 = v_new
                .iter()
                .zip(v.iter())
                .map(|(a, b)| (a - b).powi(2))
                .sum::<f64>()
                .sqrt();
            v = v_new;
            if diff < 1e-12 {
                break;
            }
        }
        for i in 0..n {
            result[i][r] = v[i];
        }
        found.push(v);
    }
    result
}
/// Higher-order SVD (HOSVD) Tucker decomposition of a 3-way tensor.
///
/// Computes a Tucker-(r0, r1, r2) approximation using sequentially truncated
/// SVD on each unfolding.
///
/// # Arguments
/// * `tensor` – 3-way `DenseTensor`
/// * `ranks`  – target rank for each mode (r0, r1, r2)
#[allow(dead_code)]
pub fn tucker_hosvd(tensor: &DenseTensor, ranks: [usize; 3]) -> TuckerDecomposition {
    assert_eq!(
        tensor.shape.len(),
        3,
        "tucker_hosvd requires a rank-3 tensor"
    );
    let [r0, r1, r2] = ranks;
    let x0 = tensor.mode_n_unfold(0);
    let x1 = tensor.mode_n_unfold(1);
    let x2 = tensor.mode_n_unfold(2);
    let u0 = truncated_svd_left(&x0, r0);
    let u1 = truncated_svd_left(&x1, r1);
    let u2 = truncated_svd_left(&x2, r2);
    let n0 = tensor.shape[0];
    let n1 = tensor.shape[1];
    let n2 = tensor.shape[2];
    let mut core = DenseTensor::zeros(&[r0, r1, r2]);
    for i in 0..r0 {
        for j in 0..r1 {
            for k in 0..r2 {
                let mut val = 0.0;
                for ii in 0..n0 {
                    for jj in 0..n1 {
                        for kk in 0..n2 {
                            val += u0[ii][i] * u1[jj][j] * u2[kk][k] * tensor.get(&[ii, jj, kk]);
                        }
                    }
                }
                core.set(&[i, j, k], val);
            }
        }
    }
    TuckerDecomposition { core, u0, u1, u2 }
}
/// Reconstruct a 3-way tensor from a Tucker decomposition.
#[allow(dead_code)]
pub fn tucker_reconstruct(td: &TuckerDecomposition) -> DenseTensor {
    let r0 = td.core.shape[0];
    let r1 = td.core.shape[1];
    let r2 = td.core.shape[2];
    let n0 = td.u0.len();
    let n1 = td.u1.len();
    let n2 = td.u2.len();
    let mut out = DenseTensor::zeros(&[n0, n1, n2]);
    for i in 0..n0 {
        for j in 0..n1 {
            for k in 0..n2 {
                let mut val = 0.0;
                for ii in 0..r0 {
                    for jj in 0..r1 {
                        for kk in 0..r2 {
                            val += td.u0[i][ii]
                                * td.u1[j][jj]
                                * td.u2[k][kk]
                                * td.core.get(&[ii, jj, kk]);
                        }
                    }
                }
                out.set(&[i, j, k], val);
            }
        }
    }
    out
}
/// Compute the bilinear form v^T A w for a Tensor2.
///
/// Returns the scalar `Σ_ij v_i A_ij w_j`.
#[allow(dead_code)]
pub fn bilinear_form(a: &Tensor2, v: &[f64; 3], w: &[f64; 3]) -> f64 {
    let mut acc = 0.0f64;
    for i in 0..3 {
        for j in 0..3 {
            acc += v[i] * a.data[i][j] * w[j];
        }
    }
    acc
}
/// Compute the quadratic form v^T A v for a Tensor2.
#[allow(dead_code)]
pub fn quadratic_form(a: &Tensor2, v: &[f64; 3]) -> f64 {
    bilinear_form(a, v, v)
}
/// Outer product of two 3-vectors: result_ij = a_i b_j (returns Tensor2).
#[allow(dead_code)]
pub fn vec_outer(a: &[f64; 3], b: &[f64; 3]) -> Tensor2 {
    let mut d = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            d[i][j] = a[i] * b[j];
        }
    }
    Tensor2 { data: d }
}
/// Tensor contraction: contract modes `mode_a` of `a` and `mode_b` of `b`
/// over their shared dimension.
///
/// For rank-2 tensors this reduces to matrix–matrix multiplication when
/// `mode_a == 1` (col) and `mode_b == 0` (row).
#[allow(dead_code)]
pub fn tensor_contraction(
    a: &DenseTensor,
    mode_a: usize,
    b: &DenseTensor,
    mode_b: usize,
) -> DenseTensor {
    let dim_a = a.shape[mode_a];
    let dim_b = b.shape[mode_b];
    assert_eq!(dim_a, dim_b, "contracted dimensions must match");
    let mut out_shape: Vec<usize> = Vec::new();
    for (k, &s) in a.shape.iter().enumerate() {
        if k != mode_a {
            out_shape.push(s);
        }
    }
    for (k, &s) in b.shape.iter().enumerate() {
        if k != mode_b {
            out_shape.push(s);
        }
    }
    if out_shape.is_empty() {
        out_shape.push(1);
    }
    let n_out: usize = out_shape.iter().product();
    let mut out_data = vec![0.0f64; n_out];
    let strides_a = a.strides();
    let strides_b = b.strides();
    let out_rank = out_shape.len();
    let mut out_strides = vec![1usize; out_rank];
    for k in (0..out_rank.saturating_sub(1)).rev() {
        out_strides[k] = out_strides[k + 1] * out_shape[k + 1];
    }
    let total_a: usize = a.data.len();
    let total_b: usize = b.data.len();
    for fa in 0..total_a {
        let mut tmp = fa;
        let mut midx_a = vec![0usize; a.shape.len()];
        for k in 0..a.shape.len() {
            midx_a[k] = tmp / strides_a[k];
            tmp %= strides_a[k];
        }
        let contract_val_a = midx_a[mode_a];
        for fb in 0..total_b {
            let mut tmp2 = fb;
            let mut midx_b = vec![0usize; b.shape.len()];
            for k in 0..b.shape.len() {
                midx_b[k] = tmp2 / strides_b[k];
                tmp2 %= strides_b[k];
            }
            if midx_b[mode_b] != contract_val_a {
                continue;
            }
            let mut oidx = vec![0usize; out_rank];
            let mut oi = 0;
            for k in 0..a.shape.len() {
                if k != mode_a {
                    oidx[oi] = midx_a[k];
                    oi += 1;
                }
            }
            for k in 0..b.shape.len() {
                if k != mode_b {
                    oidx[oi] = midx_b[k];
                    oi += 1;
                }
            }
            let fo: usize = oidx
                .iter()
                .zip(out_strides.iter())
                .map(|(&i, &s)| i * s)
                .sum();
            out_data[fo] += a.data[fa] * b.data[fb];
        }
    }
    DenseTensor {
        shape: out_shape,
        data: out_data,
    }
}
/// Compute the general outer product of two DenseTensors: result has rank = rank(a) + rank(b).
#[allow(dead_code)]
pub fn tensor_outer(a: &DenseTensor, b: &DenseTensor) -> DenseTensor {
    let mut shape = a.shape.clone();
    shape.extend_from_slice(&b.shape);
    let n_out: usize = shape.iter().product();
    let mut data = vec![0.0f64; n_out];
    let stride_b: usize = b.data.len();
    for (ia, &va) in a.data.iter().enumerate() {
        for (ib, &vb) in b.data.iter().enumerate() {
            data[ia * stride_b + ib] = va * vb;
        }
    }
    DenseTensor { shape, data }
}
/// Check whether a DenseTensor is symmetric under swap of two given modes.
#[allow(dead_code)]
pub fn is_symmetric_modes(t: &DenseTensor, mode1: usize, mode2: usize, tol: f64) -> bool {
    assert_eq!(
        t.shape[mode1], t.shape[mode2],
        "modes must have equal dimension"
    );
    let n = t.shape[mode1];
    let strides = t.strides();
    let total = t.data.len();
    for flat in 0..total {
        let mut tmp = flat;
        let mut midx = vec![0usize; t.shape.len()];
        for k in 0..t.shape.len() {
            midx[k] = tmp / strides[k];
            tmp %= strides[k];
        }
        if midx[mode1] >= midx[mode2] {
            continue;
        }
        let mut midx2 = midx.clone();
        midx2[mode1] = midx[mode2];
        midx2[mode2] = midx[mode1];
        let _n = n;
        let flat2: usize = midx2.iter().zip(strides.iter()).map(|(&i, &s)| i * s).sum();
        if (t.data[flat] - t.data[flat2]).abs() > tol {
            return false;
        }
    }
    true
}
/// Build a fully symmetric DenseTensor by averaging over all permutations of indices
/// (only for rank-2 and rank-3 tensors for now; rank-2 reduces to symmetrising a matrix).
#[allow(dead_code)]
pub fn symmetrize_tensor(t: &DenseTensor) -> DenseTensor {
    let rank = t.shape.len();
    match rank {
        2 => {
            let n = t.shape[0];
            assert_eq!(t.shape[1], n, "rank-2 symmetrize requires square tensor");
            let mut out = DenseTensor::zeros(&t.shape);
            for i in 0..n {
                for j in 0..n {
                    let v = 0.5 * (t.get(&[i, j]) + t.get(&[j, i]));
                    out.set(&[i, j], v);
                }
            }
            out
        }
        3 => {
            let n = t.shape[0];
            assert!(
                t.shape.iter().all(|&s| s == n),
                "rank-3 symmetrize requires cubic tensor"
            );
            let mut out = DenseTensor::zeros(&t.shape);
            let perms: [[usize; 3]; 6] = [
                [0, 1, 2],
                [0, 2, 1],
                [1, 0, 2],
                [1, 2, 0],
                [2, 0, 1],
                [2, 1, 0],
            ];
            for i in 0..n {
                for j in 0..n {
                    for k in 0..n {
                        let idx_orig = [i, j, k];
                        let s: f64 = perms
                            .iter()
                            .map(|p| {
                                let perm_idx = [idx_orig[p[0]], idx_orig[p[1]], idx_orig[p[2]]];
                                t.get(&perm_idx)
                            })
                            .sum::<f64>()
                            / 6.0;
                        out.set(&[i, j, k], s);
                    }
                }
            }
            out
        }
        _ => t.clone(),
    }
}
/// The 4th-order symmetric identity tensor I^S_ijkl = (delta_ik delta_jl + delta_il delta_jk)/2.
///
/// This is the projector onto symmetric second-order tensors:
/// I^S : A = sym(A) = (A + A^T)/2.
#[allow(dead_code)]
pub fn fourth_order_symmetric_identity() -> Tensor4 {
    let mut t = Tensor4::zero();
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                for l in 0..3 {
                    let v = 0.5
                        * (if i == k && j == l { 1.0 } else { 0.0 }
                            + if i == l && j == k { 1.0 } else { 0.0 });
                    t.data[i][j][k][l] = v;
                }
            }
        }
    }
    t
}
/// The 4th-order skew (antisymmetric) identity tensor:
/// I^A_ijkl = (delta_ik delta_jl - delta_il delta_jk)/2.
#[allow(dead_code)]
pub fn fourth_order_skew_identity() -> Tensor4 {
    let mut t = Tensor4::zero();
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                for l in 0..3 {
                    let v = 0.5
                        * (if i == k && j == l { 1.0 } else { 0.0 }
                            - if i == l && j == k { 1.0 } else { 0.0 });
                    t.data[i][j][k][l] = v;
                }
            }
        }
    }
    t
}
/// Compute the double inner product C :: T for a 4th-order elasticity tensor C
/// and a symmetric second-order tensor T.  Returns sigma_ij = C_ijkl T_kl.
#[allow(dead_code)]
pub fn elasticity_stress(c: &Tensor4, strain: &Tensor2) -> Tensor2 {
    c.double_contract_2(strain)
}
/// Compliance tensor from stiffness: S = C^{-1} in Voigt form (6×6 matrix inversion).
///
/// Returns `None` if the Voigt matrix is singular (det ≈ 0).
#[allow(dead_code)]
pub fn compliance_from_stiffness(c: &Tensor4) -> Option<[[f64; 6]; 6]> {
    let m = KelvinTensor::from_tensor4(c);
    invert_6x6(&m)
}
/// Invert a 6×6 matrix using Gauss–Jordan elimination.
#[allow(dead_code)]
pub fn invert_6x6(m: &[[f64; 6]; 6]) -> Option<[[f64; 6]; 6]> {
    let mut a = *m;
    let mut inv = [[0.0f64; 6]; 6];
    for i in 0..6 {
        inv[i][i] = 1.0;
    }
    for col in 0..6 {
        let mut pivot_row = col;
        let mut max_val = a[col][col].abs();
        for row in (col + 1)..6 {
            if a[row][col].abs() > max_val {
                max_val = a[row][col].abs();
                pivot_row = row;
            }
        }
        if max_val < 1e-14 {
            return None;
        }
        a.swap(col, pivot_row);
        inv.swap(col, pivot_row);
        let piv = a[col][col];
        for j in 0..6 {
            a[col][j] /= piv;
            inv[col][j] /= piv;
        }
        for row in 0..6 {
            if row == col {
                continue;
            }
            let factor = a[row][col];
            for j in 0..6 {
                let av = a[col][j];
                a[row][j] -= factor * av;
                let iv = inv[col][j];
                inv[row][j] -= factor * iv;
            }
        }
    }
    Some(inv)
}
/// Build the Eshelby inclusion tensor for an isotropic matrix.
///
/// Valid for a spherical inclusion (aspect ratio = 1).
/// Returns the Eshelby tensor S_ijkl.
///
/// Reference: J.D. Eshelby (1957), nu = Poisson's ratio.
#[allow(dead_code)]
pub fn eshelby_sphere(nu: f64) -> Tensor4 {
    let mut s = Tensor4::zero();
    let s_iiii = (7.0 - 5.0 * nu) / (15.0 * (1.0 - nu));
    let s_iijj = (5.0 * nu - 1.0) / (15.0 * (1.0 - nu));
    let s_ijij = (4.0 - 5.0 * nu) / (15.0 * (1.0 - nu));
    for i in 0..3 {
        s.data[i][i][i][i] = s_iiii;
    }
    for i in 0..3 {
        for j in 0..3 {
            if i != j {
                s.data[i][i][j][j] = s_iijj;
                s.data[i][j][i][j] = s_ijij;
                s.data[i][j][j][i] = s_ijij;
            }
        }
    }
    s
}
#[cfg(test)]
mod tests {
    use super::*;
    const EPS: f64 = 1e-12;
    fn approx_eq(a: f64, b: f64) -> bool {
        (a - b).abs() < EPS
    }
    fn tensor_approx_eq(a: &Tensor2, b: &Tensor2) -> bool {
        for i in 0..3 {
            for j in 0..3 {
                if !approx_eq(a.data[i][j], b.data[i][j]) {
                    return false;
                }
            }
        }
        true
    }
    #[test]
    fn test_identity_double_contract_equals_trace() {
        let a = Tensor2::new([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]);
        let id = Tensor2::identity();
        let dc = id.double_contract(&a);
        assert!(
            approx_eq(dc, a.trace()),
            "I:A should equal tr(A), got {dc} vs {}",
            a.trace()
        );
    }
    #[test]
    fn test_det_identity_is_one() {
        let id = Tensor2::identity();
        assert!(approx_eq(id.det(), 1.0), "det(I) should be 1");
    }
    #[test]
    fn test_inverse_of_identity_is_identity() {
        let id = Tensor2::identity();
        let inv = id.inverse().expect("identity is invertible");
        assert!(
            tensor_approx_eq(&inv, &Tensor2::identity()),
            "inv(I) should be I"
        );
    }
    #[test]
    fn test_outer_product_e1_e2() {
        let t = Tensor2::outer_product([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let mut expected = [[0.0f64; 3]; 3];
        expected[0][1] = 1.0;
        assert!(tensor_approx_eq(&t, &Tensor2::new(expected)));
    }
    #[test]
    fn test_von_mises_pure_shear() {
        let tau = 1.0f64;
        let t = Tensor2::new([[0.0, tau, 0.0], [tau, 0.0, 0.0], [0.0, 0.0, 0.0]]);
        let vm = t.von_mises();
        let expected = 3.0f64.sqrt() * tau;
        assert!(
            (vm - expected).abs() < 1e-10,
            "von Mises for pure shear: got {vm}, expected {expected}"
        );
    }
    #[test]
    fn test_isotropic_c_double_contract_identity() {
        let lambda = 1.0f64;
        let mu = 1.0f64;
        let c = Tensor4::isotropic(lambda, mu);
        let id = Tensor2::identity();
        let result = c.double_contract_2(&id);
        let expected = Tensor2::identity().scale(3.0 * lambda + 2.0 * mu);
        assert!(
            tensor_approx_eq(&result, &expected),
            "C(λ=1,μ=1):I should be (3λ+2μ)I = 5I"
        );
    }
    #[test]
    fn test_eigenvalues_diagonal_matrix() {
        let t = Tensor2::new([[2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 5.0]]);
        let ev = t.eigenvalues_symmetric();
        assert!(
            (ev[0] - 2.0).abs() < 1e-10
                && (ev[1] - 3.0).abs() < 1e-10
                && (ev[2] - 5.0).abs() < 1e-10,
            "eigenvalues should be [2,3,5], got {ev:?}"
        );
    }
    #[test]
    fn test_deviatoric_trace_is_zero() {
        let t = Tensor2::new([[3.0, 1.0, 0.5], [1.0, 2.0, 0.3], [0.5, 0.3, 4.0]]);
        let dev = t.deviatoric();
        assert!(
            dev.trace().abs() < EPS,
            "trace of deviatoric should be 0, got {}",
            dev.trace()
        );
    }
    #[test]
    fn test_voigt_round_trip() {
        let t = Tensor2::new([[1.0, 2.0, 3.0], [2.0, 4.0, 5.0], [3.0, 5.0, 6.0]]);
        let v = VoigtTensor::from_tensor2(&t);
        let t2 = VoigtTensor::to_tensor2(&v);
        assert!(tensor_approx_eq(&t, &t2), "Voigt round-trip failed");
    }
    #[test]
    fn test_contract_vec() {
        let id = Tensor2::identity();
        let v = [1.0, 2.0, 3.0];
        let result = id.contract_vec(v);
        assert!(approx_eq(result[0], 1.0));
        assert!(approx_eq(result[1], 2.0));
        assert!(approx_eq(result[2], 3.0));
    }
    #[test]
    fn test_dyadic_product() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        let t = Tensor2::dyadic(a, b);
        assert!(approx_eq(t.data[0][1], 1.0));
        assert!(approx_eq(t.data[0][0], 0.0));
    }
    #[test]
    fn test_rotate_identity() {
        let a = Tensor2::new([[1.0, 2.0, 0.0], [2.0, 3.0, 0.0], [0.0, 0.0, 4.0]]);
        let r = Tensor2::identity();
        let rotated = a.rotate(&r);
        assert!(
            tensor_approx_eq(&a, &rotated),
            "Rotating by I should be identity"
        );
    }
    #[test]
    fn test_rotate_90_degrees_z() {
        let r = Tensor2::new([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]);
        let a = Tensor2::diagonal([1.0, 0.0, 0.0]);
        let rotated = a.rotate(&r);
        assert!(
            approx_eq(rotated.data[1][1], 1.0),
            "After 90° rotation, sigma_xx should become sigma_yy"
        );
        assert!(
            rotated.data[0][0].abs() < EPS,
            "sigma_xx should be 0 after rotation"
        );
    }
    #[test]
    fn test_invariant_i1() {
        let t = Tensor2::diagonal([1.0, 2.0, 3.0]);
        assert!(approx_eq(t.invariant_i1(), 6.0));
    }
    #[test]
    fn test_invariant_i2_diagonal() {
        let t = Tensor2::diagonal([1.0, 2.0, 3.0]);
        assert!(
            approx_eq(t.invariant_i2(), 11.0),
            "I2 = {}",
            t.invariant_i2()
        );
    }
    #[test]
    fn test_invariant_i3_diagonal() {
        let t = Tensor2::diagonal([1.0, 2.0, 3.0]);
        assert!(approx_eq(t.invariant_i3(), 6.0), "I3 = det = 1*2*3 = 6");
    }
    #[test]
    fn test_principal_invariants() {
        let t = Tensor2::identity();
        let inv = t.principal_invariants();
        assert!(approx_eq(inv[0], 3.0));
        assert!(approx_eq(inv[1], 3.0));
        assert!(approx_eq(inv[2], 1.0));
    }
    #[test]
    fn test_decompose_dev_hydro() {
        let t = Tensor2::diagonal([3.0, 6.0, 9.0]);
        let (dev, h) = t.decompose_dev_hydro();
        assert!(approx_eq(h, 6.0), "hydrostatic = {h}");
        assert!(
            dev.trace().abs() < EPS,
            "deviatoric trace = {}",
            dev.trace()
        );
    }
    #[test]
    fn test_from_voigt() {
        let v = [1.0, 2.0, 3.0, 0.5, 0.3, 0.1];
        let t = Tensor2::from_voigt(v);
        assert!(approx_eq(t.data[0][0], 1.0));
        assert!(approx_eq(t.data[1][1], 2.0));
        assert!(approx_eq(t.data[0][1], 0.5));
        assert!(approx_eq(t.data[1][0], 0.5));
    }
    #[test]
    fn test_to_voigt() {
        let t = Tensor2::new([[1.0, 0.5, 0.1], [0.5, 2.0, 0.3], [0.1, 0.3, 3.0]]);
        let v = t.to_voigt();
        assert!(approx_eq(v[0], 1.0));
        assert!(approx_eq(v[3], 0.5));
    }
    #[test]
    fn test_effective_strain() {
        let t = Tensor2::zero();
        assert!(
            t.effective_strain() < EPS,
            "zero tensor should give zero strain"
        );
    }
    #[test]
    fn test_diagonal_constructor() {
        let t = Tensor2::diagonal([1.0, 2.0, 3.0]);
        assert!(approx_eq(t.data[0][0], 1.0));
        assert!(approx_eq(t.data[1][1], 2.0));
        assert!(approx_eq(t.data[2][2], 3.0));
        assert!(approx_eq(t.data[0][1], 0.0));
    }
    #[test]
    fn test_is_symmetric() {
        let t = Tensor2::new([[1.0, 2.0, 3.0], [2.0, 4.0, 5.0], [3.0, 5.0, 6.0]]);
        assert!(t.is_symmetric(1e-14));
        let t2 = Tensor2::new([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]);
        assert!(!t2.is_symmetric(1e-14));
    }
    #[test]
    fn test_matrix_exp_identity() {
        let zero = Tensor2::zero();
        let result = zero.matrix_exp_approx();
        assert!(
            tensor_approx_eq(&result, &Tensor2::identity()),
            "exp(0) = I"
        );
    }
    #[test]
    fn test_deviatoric_projector() {
        let p = Tensor4::deviatoric_projector();
        let id = Tensor2::identity();
        let result = p.double_contract_2(&id);
        assert!(
            result.norm() < 1e-10,
            "P_dev : I should be zero: norm = {}",
            result.norm()
        );
    }
    #[test]
    fn test_tensor4_to_voigt_matrix_isotropic() {
        let lambda = 1.0;
        let mu = 1.0;
        let c = Tensor4::isotropic(lambda, mu);
        let m = c.to_voigt_matrix();
        assert!(approx_eq(m[0][0], 3.0), "C_1111 = {}", m[0][0]);
        assert!(approx_eq(m[0][1], 1.0), "C_1122 = {}", m[0][1]);
        assert!(approx_eq(m[3][3], 1.0), "C_1212 = {}", m[3][3]);
    }
    #[test]
    fn test_tensor4_scale() {
        let c = Tensor4::isotropic(1.0, 1.0);
        let c2 = c.scale(2.0);
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        assert!(approx_eq(c2.data[i][j][k][l], 2.0 * c.data[i][j][k][l]));
                    }
                }
            }
        }
    }
    #[test]
    fn test_tensor4_add() {
        let c1 = Tensor4::isotropic(1.0, 1.0);
        let c2 = Tensor4::isotropic(2.0, 3.0);
        let c3 = c1.add(&c2);
        let expected = Tensor4::isotropic(3.0, 4.0);
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        assert!(approx_eq(c3.data[i][j][k][l], expected.data[i][j][k][l]));
                    }
                }
            }
        }
    }
    #[test]
    fn test_engineering_strain_round_trip() {
        let v = [0.01, 0.02, 0.03, 0.005, 0.003, 0.001];
        let eng = VoigtTensor::to_engineering_strain(&v);
        let back = VoigtTensor::from_engineering_strain(&eng);
        for i in 0..6 {
            assert!(approx_eq(v[i], back[i]), "round trip failed at {i}");
        }
    }
    #[test]
    fn test_von_mises_voigt() {
        let v = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        let vm = VoigtTensor::von_mises_voigt(&v);
        let expected = 3.0_f64.sqrt();
        assert!(
            (vm - expected).abs() < 1e-10,
            "von Mises = {vm}, expected {expected}"
        );
    }
    #[test]
    fn test_hydrostatic_pressure() {
        let v = [-100.0, -100.0, -100.0, 0.0, 0.0, 0.0];
        let p = VoigtTensor::hydrostatic_pressure(&v);
        assert!(approx_eq(p, 100.0), "pressure = {p}, expected 100.0");
    }
    #[test]
    fn test_sym_and_skew_add_to_original() {
        let a = Tensor2::new([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]);
        let sym = a.sym();
        let skew = a.skew();
        let reconstructed = sym.add(&skew);
        assert!(
            tensor_approx_eq(&a, &reconstructed),
            "sym + skew should equal original"
        );
    }
    #[test]
    fn test_tensor4_has_minor_symmetry_isotropic() {
        let c = Tensor4::isotropic(1.0, 1.0);
        assert!(
            c.has_minor_symmetry(1e-12),
            "isotropic tensor should have minor symmetry"
        );
    }
    #[test]
    fn test_tensor4_has_major_symmetry_isotropic() {
        let c = Tensor4::isotropic(2.0, 3.0);
        assert!(
            c.has_major_symmetry(1e-12),
            "isotropic tensor should have major symmetry"
        );
    }
    #[test]
    fn test_tensor4_symmetrize_minor_idempotent() {
        let c = Tensor4::isotropic(1.0, 1.0);
        let cs = c.symmetrize_minor();
        let cs2 = cs.symmetrize_minor();
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        assert!(approx_eq(cs.data[i][j][k][l], cs2.data[i][j][k][l]));
                    }
                }
            }
        }
    }
    #[test]
    fn test_tensor4_symmetrize_major_idempotent() {
        let c = Tensor4::isotropic(1.0, 1.0);
        let cs = c.symmetrize_major();
        let cs2 = cs.symmetrize_major();
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        assert!(approx_eq(cs.data[i][j][k][l], cs2.data[i][j][k][l]));
                    }
                }
            }
        }
    }
    #[test]
    fn test_tensor4_double_contract_4_identity() {
        let c = Tensor4::isotropic(1.0, 1.0);
        let id4 = Tensor4::identity();
        let result = c.double_contract_4(&id4);
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        assert!(
                            approx_eq(result.data[i][j][k][l], c.data[i][j][k][l]),
                            "C::I4 should equal C at [{i}][{j}][{k}][{l}]"
                        );
                    }
                }
            }
        }
    }
    #[test]
    fn test_tensor4_rotate_identity_unchanged() {
        let c = Tensor4::isotropic(1.0, 1.0);
        let r = Tensor2::identity();
        let c_rot = c.rotate(&r);
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        assert!(
                            (c_rot.data[i][j][k][l] - c.data[i][j][k][l]).abs() < 1e-10,
                            "rotating by I should be unchanged"
                        );
                    }
                }
            }
        }
    }
    #[test]
    fn test_tensor4_norm_positive() {
        let c = Tensor4::isotropic(1.0, 1.0);
        assert!(c.norm() > 0.0, "Frobenius norm should be positive");
    }
    #[test]
    fn test_tensor4_sub() {
        let c1 = Tensor4::isotropic(2.0, 2.0);
        let c2 = Tensor4::isotropic(1.0, 1.0);
        let diff = c1.sub(&c2);
        let expected = Tensor4::isotropic(1.0, 1.0);
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        assert!(approx_eq(diff.data[i][j][k][l], expected.data[i][j][k][l]));
                    }
                }
            }
        }
    }
    #[test]
    fn test_tensor4_deviatoric_plus_volumetric_is_identity_sym() {
        let pd = Tensor4::deviatoric_projector();
        let pv = Tensor4::volumetric_projector();
        let sum = pd.add(&pv);
        let id_sym = Tensor4::identity_sym();
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        assert!(
                            (sum.data[i][j][k][l] - id_sym.data[i][j][k][l]).abs() < 1e-12,
                            "P_dev + P_vol != I_sym at [{i}][{j}][{k}][{l}]"
                        );
                    }
                }
            }
        }
    }
    #[test]
    fn test_kelvin_from_tensor4_isotropic_diagonal() {
        let c = Tensor4::isotropic(1.0, 1.0);
        let km = KelvinTensor::from_tensor4(&c);
        assert!(
            (km[0][0] - 3.0).abs() < 1e-10,
            "K[0][0] = {} expected 3",
            km[0][0]
        );
        assert!(
            (km[3][3] - 2.0).abs() < 1e-10,
            "K[3][3] = {} expected 2",
            km[3][3]
        );
    }
    #[test]
    fn test_kelvin_round_trip() {
        let c = Tensor4::isotropic(2.0, 3.0);
        let km = KelvinTensor::from_tensor4(&c);
        let c2 = KelvinTensor::to_tensor4(&km);
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        assert!(
                            (c.data[i][j][k][l] - c2.data[i][j][k][l]).abs() < 1e-10,
                            "Kelvin round-trip failed at [{i}][{j}][{k}][{l}]"
                        );
                    }
                }
            }
        }
    }
    #[test]
    fn test_kelvin_stress_round_trip() {
        let sigma = Tensor2::new([[1.0, 0.5, 0.3], [0.5, 2.0, 0.2], [0.3, 0.2, 3.0]]);
        let kv = KelvinTensor::stress_to_kelvin(&sigma);
        let sigma2 = KelvinTensor::kelvin_to_stress(&kv);
        assert!(
            tensor_approx_eq(&sigma, &sigma2),
            "Kelvin stress round-trip failed"
        );
    }
    #[test]
    fn test_kelvin_matvec_isotropic() {
        let lambda = 1.0_f64;
        let mu = 1.0_f64;
        let c = Tensor4::isotropic(lambda, mu);
        let km = KelvinTensor::from_tensor4(&c);
        let id_kelvin = KelvinTensor::stress_to_kelvin(&Tensor2::identity());
        let result_kelvin = KelvinTensor::matvec(&km, &id_kelvin);
        let result = KelvinTensor::kelvin_to_stress(&result_kelvin);
        let expected_val = 3.0 * lambda + 2.0 * mu;
        assert!((result.data[0][0] - expected_val).abs() < 1e-10);
        assert!((result.data[1][1] - expected_val).abs() < 1e-10);
        assert!((result.data[0][1]).abs() < 1e-10);
    }
    #[test]
    fn test_rotation_z_is_orthogonal() {
        use std::f64::consts::PI;
        let r = TensorBasis::rotation_z(PI / 6.0);
        let rt = r.transpose();
        let rrt = r.dot(&rt);
        assert!(
            tensor_approx_eq(&rrt, &Tensor2::identity()),
            "R*R^T should be I"
        );
    }
    #[test]
    fn test_rotation_x_is_orthogonal() {
        let r = TensorBasis::rotation_x(0.7);
        let rrt = r.dot(&r.transpose());
        assert!(tensor_approx_eq(&rrt, &Tensor2::identity()));
    }
    #[test]
    fn test_rotate_voigt_stiffness_identity_rotation() {
        let c = Tensor4::isotropic(1.0, 1.0);
        let m = c.to_voigt_matrix();
        let r = Tensor2::identity();
        let m_rot = TensorBasis::rotate_voigt_stiffness(&m, &r);
        for p in 0..6 {
            for q in 0..6 {
                assert!(
                    (m[p][q] - m_rot[p][q]).abs() < 1e-10,
                    "rotation by I should not change stiffness matrix at [{p}][{q}]"
                );
            }
        }
    }
    #[test]
    fn test_tensor4_single_contract_right_identity() {
        let c = Tensor4::isotropic(1.0, 1.0);
        let id = Tensor2::identity();
        let result = c.single_contract_right(&id);
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        assert!(
                            (result.data[i][j][k][l] - c.data[i][j][k][l]).abs() < 1e-10,
                            "C·I should equal C at [{i}][{j}][{k}][{l}]"
                        );
                    }
                }
            }
        }
    }
    #[test]
    fn test_tensor3_zero() {
        let t = Tensor3::zero();
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    assert_eq!(t.data[i][j][k], 0.0, "Tensor3::zero should be all zeros");
                }
            }
        }
    }
    #[test]
    fn test_tensor3_set_get() {
        let mut t = Tensor3::zero();
        t.data[1][2][0] = 3.125;
        assert!((t.data[1][2][0] - 3.125).abs() < 1e-12);
    }
    #[test]
    fn test_tensor3_levi_civita() {
        let eps = Tensor3::levi_civita();
        assert!((eps.data[0][1][2] - 1.0).abs() < EPS, "ε_012 should be 1");
        assert!((eps.data[0][2][1] + 1.0).abs() < EPS, "ε_021 should be -1");
        assert!((eps.data[1][0][2] + 1.0).abs() < EPS, "ε_102 should be -1");
        assert!((eps.data[1][2][0] - 1.0).abs() < EPS, "ε_120 should be 1");
        assert!((eps.data[2][0][1] - 1.0).abs() < EPS, "ε_201 should be 1");
        assert!((eps.data[2][1][0] + 1.0).abs() < EPS, "ε_210 should be -1");
        assert!((eps.data[0][0][0]).abs() < EPS, "ε_000 should be 0");
    }
    #[test]
    fn test_tensor3_contract_vec_with_levi_civita() {
        let eps = Tensor3::levi_civita();
        let v = [1.0_f64, 0.0, 0.0];
        let result = eps.contract_last(&v);
        assert!(
            (result.data[1][2] - 1.0).abs() < EPS,
            "ε_12k*e_0_k = ε_120 = 1"
        );
        assert!(
            (result.data[2][1] + 1.0).abs() < EPS,
            "ε_21k*e_0_k = ε_210 = -1"
        );
    }
    #[test]
    fn test_tensor3_scale() {
        let eps = Tensor3::levi_civita();
        let scaled = eps.scale(2.0);
        assert!(
            (scaled.data[0][1][2] - 2.0).abs() < EPS,
            "scaled ε_012 should be 2"
        );
    }
    #[test]
    fn test_tensor3_add_zero() {
        let eps = Tensor3::levi_civita();
        let zero = Tensor3::zero();
        let sum = eps.add(&zero);
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    assert!(
                        (sum.data[i][j][k] - eps.data[i][j][k]).abs() < EPS,
                        "add with zero should not change tensor at [{i}][{j}][{k}]"
                    );
                }
            }
        }
    }
    #[test]
    fn test_tensor3_totally_antisymmetric_levi_civita() {
        let eps = Tensor3::levi_civita();
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    let e_ijk = eps.data[i][j][k];
                    let e_jki = eps.data[j][k][i];
                    let e_kij = eps.data[k][i][j];
                    assert!(
                        (e_ijk - e_jki).abs() < EPS,
                        "ε_ijk != ε_jki for ({i},{j},{k})"
                    );
                    assert!(
                        (e_ijk - e_kij).abs() < EPS,
                        "ε_ijk != ε_kij for ({i},{j},{k})"
                    );
                }
            }
        }
    }
    #[test]
    fn test_mandel_normal_components_unchanged() {
        let sigma = Tensor2::new([[1.0, 0.5, 0.3], [0.5, 2.0, 0.2], [0.3, 0.2, 3.0]]);
        let mv = MandelNotation::from_tensor2(&sigma);
        assert!((mv[0] - 1.0).abs() < EPS, "Mandel[0] = σ_11 = 1");
        assert!((mv[1] - 2.0).abs() < EPS, "Mandel[1] = σ_22 = 2");
        assert!((mv[2] - 3.0).abs() < EPS, "Mandel[2] = σ_33 = 3");
    }
    #[test]
    fn test_mandel_shear_scaled_by_sqrt2() {
        let sigma = Tensor2::new([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]]);
        let mv = MandelNotation::from_tensor2(&sigma);
        let s2 = std::f64::consts::SQRT_2;
        assert!((mv[3] - s2).abs() < EPS, "Mandel[3] = √2 * σ_12");
    }
    #[test]
    fn test_mandel_round_trip() {
        let sigma = Tensor2::new([[1.0, 0.5, 0.3], [0.5, 2.0, 0.2], [0.3, 0.2, 3.0]]);
        let mv = MandelNotation::from_tensor2(&sigma);
        let sigma2 = MandelNotation::to_tensor2(&mv);
        assert!(
            tensor_approx_eq(&sigma, &sigma2),
            "Mandel round-trip failed"
        );
    }
    #[test]
    fn test_mandel_inner_product_preserved() {
        let a = Tensor2::new([[1.0, 0.5, 0.0], [0.5, 2.0, 0.0], [0.0, 0.0, 3.0]]);
        let b = Tensor2::new([[2.0, 0.1, 0.0], [0.1, 1.0, 0.0], [0.0, 0.0, 0.5]]);
        let ab_direct = a.double_contract(&b);
        let ma = MandelNotation::from_tensor2(&a);
        let mb = MandelNotation::from_tensor2(&b);
        let ab_mandel: f64 = ma.iter().zip(mb.iter()).map(|(x, y)| x * y).sum();
        assert!(
            (ab_direct - ab_mandel).abs() < 1e-10,
            "Mandel inner product should match direct double contraction: {} vs {}",
            ab_direct,
            ab_mandel
        );
    }
    #[test]
    fn test_mandel_norm_equals_frobenius() {
        let a = Tensor2::new([[1.0, 0.5, 0.3], [0.5, 2.0, 0.2], [0.3, 0.2, 3.0]]);
        let frob_sq: f64 = (0..3)
            .flat_map(|i| (0..3).map(move |j| a.data[i][j] * a.data[i][j]))
            .sum();
        let mv = MandelNotation::from_tensor2(&a);
        let mandel_sq: f64 = mv.iter().map(|x| x * x).sum();
        assert!(
            (frob_sq - mandel_sq).abs() < 1e-10,
            "Mandel norm² should equal Frobenius norm²: {} vs {}",
            frob_sq,
            mandel_sq
        );
    }
    #[test]
    fn test_spectral_decomp_identity() {
        let id = Tensor2::identity();
        let evals = id.eigenvalues_symmetric();
        for &e in &evals {
            assert!(
                (e - 1.0).abs() < 1e-8,
                "Identity eigenvalue should be 1, got {e}"
            );
        }
    }
    #[test]
    fn test_spectral_decomp_diagonal_sorted() {
        let d = Tensor2::diagonal([3.0, 1.0, 2.0]);
        let evals = d.eigenvalues_symmetric();
        assert!(
            evals[0] <= evals[1] && evals[1] <= evals[2],
            "Eigenvalues should be sorted ascending: {:?}",
            evals
        );
        assert!((evals[0] - 1.0).abs() < 1e-8);
        assert!((evals[1] - 2.0).abs() < 1e-8);
        assert!((evals[2] - 3.0).abs() < 1e-8);
    }
    #[test]
    fn test_spectral_decomp_stress_tensor() {
        let sigma = Tensor2::new([[2.0, 1.0, 0.0], [1.0, 2.0, 0.0], [0.0, 0.0, 1.0]]);
        let evals = sigma.eigenvalues_symmetric();
        let mut sorted = evals;
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!(
            (sorted[0] - 1.0).abs() < 1e-6,
            "eval[0] = {} expected 1",
            sorted[0]
        );
        assert!(
            (sorted[1] - 1.0).abs() < 1e-6,
            "eval[1] = {} expected 1",
            sorted[1]
        );
        assert!(
            (sorted[2] - 3.0).abs() < 1e-6,
            "eval[2] = {} expected 3",
            sorted[2]
        );
    }
    #[test]
    fn test_voigt6x6_isotropic_symmetry() {
        let c = Tensor4::isotropic(1.0, 1.0);
        let m = c.to_voigt_matrix();
        for p in 0..6 {
            for q in 0..6 {
                assert!(
                    (m[p][q] - m[q][p]).abs() < 1e-10,
                    "Voigt 6x6 should be symmetric at [{p}][{q}]: {} vs {}",
                    m[p][q],
                    m[q][p]
                );
            }
        }
    }
    #[test]
    fn test_voigt6x6_positive_diagonal() {
        let c = Tensor4::isotropic(1.0, 1.0);
        let m = c.to_voigt_matrix();
        for p in 0..6 {
            assert!(
                m[p][p] > 0.0,
                "Voigt diagonal M[{p}][{p}] = {} should be > 0",
                m[p][p]
            );
        }
    }
    #[test]
    fn test_voigt6x6_isotropic_values() {
        let c = Tensor4::isotropic(1.0, 1.0);
        let m = c.to_voigt_matrix();
        assert!(
            (m[0][0] - 3.0).abs() < 1e-10,
            "C_1111 = λ+2μ = 3, got {}",
            m[0][0]
        );
        assert!(
            (m[0][1] - 1.0).abs() < 1e-10,
            "C_1122 = λ = 1, got {}",
            m[0][1]
        );
        assert!(
            (m[3][3] - 1.0).abs() < 1e-10,
            "C_1212 = μ = 1, got {}",
            m[3][3]
        );
    }
    #[test]
    fn test_double_contract_commutative_for_symmetric() {
        let a = Tensor2::new([[1.0, 0.5, 0.0], [0.5, 2.0, 0.0], [0.0, 0.0, 3.0]]);
        let b = Tensor2::new([[2.0, 0.1, 0.0], [0.1, 1.0, 0.0], [0.0, 0.0, 0.5]]);
        let ab = a.double_contract(&b);
        let ba = b.double_contract(&a);
        assert!(
            (ab - ba).abs() < 1e-10,
            "A:B should equal B:A for symmetric tensors"
        );
    }
    #[test]
    fn test_deviatoric_plus_hydrostatic_equals_original() {
        let sigma = Tensor2::new([[3.0, 1.0, 0.5], [1.0, 4.0, 0.2], [0.5, 0.2, 5.0]]);
        let (dev, p) = sigma.decompose_dev_hydro();
        let hydro = Tensor2::identity().scale(p);
        let reconstructed = dev.add(&hydro);
        assert!(
            tensor_approx_eq(&reconstructed, &sigma),
            "dev + hydro should equal original tensor"
        );
    }
    #[test]
    fn test_deviatoric_is_traceless() {
        let sigma = Tensor2::new([[3.0, 1.0, 0.5], [1.0, 4.0, 0.2], [0.5, 0.2, 5.0]]);
        let dev = sigma.deviatoric();
        let tr = dev.trace();
        assert!(
            tr.abs() < 1e-10,
            "Deviatoric tensor should be traceless, got {tr}"
        );
    }
    #[test]
    fn test_major_symmetry_isotropic() {
        let c = Tensor4::isotropic(2.0, 3.0);
        assert!(
            c.has_major_symmetry(1e-10),
            "Isotropic C should have major symmetry"
        );
    }
    #[test]
    fn test_minor_symmetry_isotropic() {
        let c = Tensor4::isotropic(2.0, 3.0);
        assert!(
            c.has_minor_symmetry(1e-10),
            "Isotropic C should have minor symmetry"
        );
    }
    #[test]
    fn test_dense_tensor_zeros() {
        let t = DenseTensor::zeros(&[2, 3, 4]);
        assert_eq!(t.shape, vec![2, 3, 4]);
        assert_eq!(t.data.len(), 24);
        assert!(t.data.iter().all(|&x| x == 0.0));
    }
    #[test]
    fn test_dense_tensor_from_data() {
        let data: Vec<f64> = (0..6).map(|x| x as f64).collect();
        let t = DenseTensor::from_data(&[2, 3], data.clone());
        assert_eq!(t.get(&[0, 0]), 0.0);
        assert_eq!(t.get(&[0, 1]), 1.0);
        assert_eq!(t.get(&[1, 2]), 5.0);
    }
    #[test]
    fn test_dense_tensor_set_get() {
        let mut t = DenseTensor::zeros(&[3, 3]);
        t.set(&[1, 2], 42.0);
        assert_eq!(t.get(&[1, 2]), 42.0);
        assert_eq!(t.get(&[0, 0]), 0.0);
    }
    #[test]
    fn test_dense_tensor_frobenius_norm() {
        let t = DenseTensor::from_data(&[2, 2], vec![1.0, 2.0, 3.0, 4.0]);
        let expected = (1.0f64 + 4.0 + 9.0 + 16.0).sqrt();
        assert!(
            (t.frobenius_norm() - expected).abs() < 1e-12,
            "Frobenius norm mismatch: {} vs {}",
            t.frobenius_norm(),
            expected
        );
    }
    #[test]
    fn test_mode_n_unfold_matrix() {
        let t = DenseTensor::from_data(&[2, 3], vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0]);
        let m0 = t.mode_n_unfold(0);
        assert_eq!(m0.len(), 2);
        assert_eq!(m0[0].len(), 3);
        let m1 = t.mode_n_unfold(1);
        assert_eq!(m1.len(), 3);
        assert_eq!(m1[0].len(), 2);
    }
    #[test]
    fn test_mode_n_fold_roundtrip() {
        let data: Vec<f64> = (0..24).map(|x| x as f64).collect();
        let t = DenseTensor::from_data(&[2, 3, 4], data.clone());
        for mode in 0..3 {
            let unfolded = t.mode_n_unfold(mode);
            let refolded = DenseTensor::mode_n_fold(&unfolded, mode, &[2, 3, 4]);
            for (a, b) in t.data.iter().zip(refolded.data.iter()) {
                assert!(
                    (a - b).abs() < 1e-12,
                    "mode_{mode} fold round-trip failed: {a} vs {b}"
                );
            }
        }
    }
    #[test]
    fn test_einsum_matmul() {
        let a = vec![vec![1.0, 2.0], vec![3.0, 4.0]];
        let b = vec![vec![5.0, 6.0], vec![7.0, 8.0]];
        let c = einsum_2d("ij,jk->ik", &a, Some(&b));
        assert!((c[0][0] - 19.0).abs() < 1e-12);
        assert!((c[0][1] - 22.0).abs() < 1e-12);
        assert!((c[1][0] - 43.0).abs() < 1e-12);
        assert!((c[1][1] - 50.0).abs() < 1e-12);
    }
    #[test]
    fn test_einsum_frobenius_inner_product() {
        let a = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let b = vec![vec![2.0, 3.0], vec![4.0, 5.0]];
        let r = einsum_2d("ij,ij->", &a, Some(&b));
        assert!((r[0][0] - 7.0).abs() < 1e-12);
    }
    #[test]
    fn test_einsum_transpose() {
        let a = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]];
        let t = einsum_2d("ij->ji", &a, None);
        assert_eq!(t.len(), 3);
        assert_eq!(t[0].len(), 2);
        assert!((t[0][0] - 1.0).abs() < 1e-12);
        assert!((t[2][1] - 6.0).abs() < 1e-12);
    }
    #[test]
    fn test_einsum_trace() {
        let a = vec![vec![1.0, 2.0], vec![3.0, 4.0]];
        let r = einsum_2d("ii->", &a, None);
        assert!(
            (r[0][0] - 5.0).abs() < 1e-12,
            "trace should be 5, got {}",
            r[0][0]
        );
    }
    #[test]
    fn test_cp_als_rank1_reconstruction() {
        let a = [1.0, 2.0];
        let b = [1.0, 2.0, 3.0];
        let c = [1.0, 2.0];
        let mut data = vec![];
        for &ai in &a {
            for &bj in &b {
                for &ck in &c {
                    data.push(ai * bj * ck);
                }
            }
        }
        let tensor = DenseTensor::from_data(&[2, 3, 2], data.clone());
        let fro_norm = tensor.frobenius_norm();
        let cp = cp_als(&tensor, 1, 500, 1e-12);
        let err = cp_reconstruction_error(&tensor, &cp.a, &cp.b, &cp.c, &cp.lambdas);
        let rel_err = err / (fro_norm + 1e-14);
        assert!(
            rel_err < 0.1,
            "CP rank-1 relative error should be < 10%, got {rel_err}"
        );
    }
    #[test]
    fn test_cp_als_returns_correct_shapes() {
        let tensor = DenseTensor::zeros(&[3, 4, 5]);
        let cp = cp_als(&tensor, 2, 10, 1e-6);
        assert_eq!(cp.a.len(), 3);
        assert_eq!(cp.b.len(), 4);
        assert_eq!(cp.c.len(), 5);
        assert_eq!(cp.lambdas.len(), 2);
    }
    #[test]
    fn test_cp_reconstruction_error_zeros() {
        let tensor = DenseTensor::zeros(&[2, 2, 2]);
        let cp = cp_als(&tensor, 1, 5, 1e-10);
        let err = cp_reconstruction_error(&tensor, &cp.a, &cp.b, &cp.c, &cp.lambdas);
        assert!(err < 1e-6, "Zero tensor CP error should be tiny, got {err}");
    }
    #[test]
    fn test_tucker_hosvd_shape() {
        let data: Vec<f64> = (0..24).map(|x| x as f64).collect();
        let tensor = DenseTensor::from_data(&[2, 3, 4], data);
        let td = tucker_hosvd(&tensor, [2, 2, 2]);
        assert_eq!(td.core.shape, vec![2, 2, 2]);
        assert_eq!(td.u0.len(), 2);
        assert_eq!(td.u0[0].len(), 2);
        assert_eq!(td.u1.len(), 3);
        assert_eq!(td.u2.len(), 4);
    }
    #[test]
    fn test_tucker_hosvd_full_rank_reconstruction() {
        let data: Vec<f64> = (0..8).map(|x| x as f64).collect();
        let tensor = DenseTensor::from_data(&[2, 2, 2], data);
        let fro = tensor.frobenius_norm();
        let td = tucker_hosvd(&tensor, [2, 2, 2]);
        let recon = tucker_reconstruct(&td);
        let err: f64 = tensor
            .data
            .iter()
            .zip(recon.data.iter())
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt();
        let rel_err = err / (fro + 1e-14);
        assert!(
            rel_err < 0.5,
            "Tucker full-rank reconstruction relative error: {rel_err} (err={err}, fro={fro})"
        );
    }
    #[test]
    fn test_tucker_reconstruct_shape() {
        let data: Vec<f64> = (0..27).map(|x| x as f64).collect();
        let tensor = DenseTensor::from_data(&[3, 3, 3], data);
        let td = tucker_hosvd(&tensor, [2, 2, 2]);
        let recon = tucker_reconstruct(&td);
        assert_eq!(recon.shape, vec![3, 3, 3]);
        assert_eq!(recon.data.len(), 27);
    }
    #[test]
    fn test_tt_core_get_set() {
        let mut core = TtCore::zeros(2, 3, 4);
        core.set(1, 2, 3, 5.0);
        assert!((core.get(1, 2, 3) - 5.0).abs() < 1e-12);
        assert!((core.get(0, 0, 0)).abs() < 1e-12);
    }
    #[test]
    fn test_tt_core_frobenius_norm() {
        let mut core = TtCore::zeros(1, 2, 1);
        core.set(0, 0, 0, 3.0);
        core.set(0, 1, 0, 4.0);
        assert!(
            (core.frobenius_norm() - 5.0).abs() < 1e-12,
            "TT core norm: expected 5, got {}",
            core.frobenius_norm()
        );
    }
    #[test]
    fn test_tt_evaluate_simple() {
        let mut c0 = TtCore::zeros(1, 2, 1);
        c0.set(0, 0, 0, 1.0);
        c0.set(0, 1, 0, 2.0);
        let mut c1 = TtCore::zeros(1, 2, 1);
        c1.set(0, 0, 0, 1.0);
        c1.set(0, 1, 0, 2.0);
        let tt = TensorTrain {
            cores: vec![c0, c1],
            shape: vec![2, 2],
        };
        assert!((tt.evaluate(&[0, 0]) - 1.0).abs() < 1e-12);
        assert!((tt.evaluate(&[0, 1]) - 2.0).abs() < 1e-12);
        assert!((tt.evaluate(&[1, 0]) - 2.0).abs() < 1e-12);
        assert!((tt.evaluate(&[1, 1]) - 4.0).abs() < 1e-12);
    }
    #[test]
    fn test_tt_from_dense_shape() {
        let data: Vec<f64> = (0..8).map(|x| x as f64).collect();
        let tensor = DenseTensor::from_data(&[2, 2, 2], data);
        let tt = TensorTrain::from_dense(&tensor, 4, 1e-10);
        assert_eq!(tt.shape, vec![2, 2, 2]);
        assert_eq!(tt.cores.len(), 3);
        assert_eq!(tt.cores[0].r_left, 1);
        assert_eq!(tt.cores.last().unwrap().r_right, 1);
    }
    #[test]
    fn test_tt_frobenius_norm_identity_like() {
        let mut c0 = TtCore::zeros(1, 2, 1);
        c0.set(0, 0, 0, 1.0);
        c0.set(0, 1, 0, 1.0);
        let mut c1 = TtCore::zeros(1, 2, 1);
        c1.set(0, 0, 0, 1.0);
        c1.set(0, 1, 0, 1.0);
        let tt = TensorTrain {
            cores: vec![c0, c1],
            shape: vec![2, 2],
        };
        let norm = tt.frobenius_norm_approx();
        assert!(
            (norm - 2.0).abs() < 1e-12,
            "TT Frobenius norm expected 2, got {norm}"
        );
    }
    #[test]
    fn test_matmul_helper() {
        let a = vec![vec![1.0, 2.0], vec![3.0, 4.0]];
        let b = vec![vec![5.0, 6.0], vec![7.0, 8.0]];
        let c = matmul(&a, &b);
        assert!((c[0][0] - 19.0).abs() < 1e-12);
        assert!((c[1][1] - 50.0).abs() < 1e-12);
    }
    #[test]
    fn test_gram_helper() {
        let a = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let g = gram(&a);
        assert!((g[0][0] - 1.0).abs() < 1e-12);
        assert!((g[0][1]).abs() < 1e-12);
        assert!((g[1][1] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_khatri_rao_shape() {
        let a = vec![vec![1.0, 2.0], vec![3.0, 4.0], vec![5.0, 6.0]];
        let b = vec![vec![7.0, 8.0], vec![9.0, 10.0]];
        let kr = khatri_rao(&a, &b);
        assert_eq!(kr.len(), 6);
        assert_eq!(kr[0].len(), 2);
    }
    #[test]
    fn test_transpose_helper() {
        let a = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]];
        let t = transpose(&a);
        assert_eq!(t.len(), 3);
        assert_eq!(t[0].len(), 2);
        assert!((t[2][1] - 6.0).abs() < 1e-12);
    }
    #[test]
    fn test_solve_ls_identity() {
        let a = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let b = vec![vec![3.0, 5.0]];
        let b_t = vec![vec![3.0], vec![5.0]];
        let x = solve_ls(&a, &b_t);
        assert!((x[0][0] - 3.0).abs() < 1e-12);
        assert!((x[1][0] - 5.0).abs() < 1e-12);
        let _ = b;
    }
}
