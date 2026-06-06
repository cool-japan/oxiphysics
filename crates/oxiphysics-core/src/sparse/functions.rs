//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{CsrMatrix, SparseLu};

/// Return the n×n identity matrix in CSR format.
pub fn identity_csr(n: usize) -> CsrMatrix {
    let rows: Vec<usize> = (0..n).collect();
    let cols: Vec<usize> = (0..n).collect();
    let vals: Vec<f64> = vec![1.0; n];
    CsrMatrix::from_triplets(n, n, &rows, &cols, &vals)
}
/// 1-D finite-difference Laplacian on `n` interior nodes.
///
/// Stencil: `(-1, 2, -1) / dx²`
pub fn laplacian_1d(n: usize, dx: f64) -> CsrMatrix {
    let inv_dx2 = 1.0 / (dx * dx);
    let mut rows = Vec::new();
    let mut cols = Vec::new();
    let mut vals = Vec::new();
    for i in 0..n {
        rows.push(i);
        cols.push(i);
        vals.push(2.0 * inv_dx2);
        if i > 0 {
            rows.push(i);
            cols.push(i - 1);
            vals.push(-inv_dx2);
        }
        if i + 1 < n {
            rows.push(i);
            cols.push(i + 1);
            vals.push(-inv_dx2);
        }
    }
    CsrMatrix::from_triplets(n, n, &rows, &cols, &vals)
}
/// 2-D finite-difference Laplacian on an `nx × ny` grid (5-point stencil).
///
/// Grid ordering: node `(i, j)` maps to global index `i * ny + j`.
pub fn laplacian_2d(nx: usize, ny: usize, dx: f64) -> CsrMatrix {
    let n = nx * ny;
    let inv_dx2 = 1.0 / (dx * dx);
    let mut rows = Vec::new();
    let mut cols = Vec::new();
    let mut vals = Vec::new();
    for i in 0..nx {
        for j in 0..ny {
            let node = i * ny + j;
            rows.push(node);
            cols.push(node);
            vals.push(4.0 * inv_dx2);
            if j > 0 {
                rows.push(node);
                cols.push(node - 1);
                vals.push(-inv_dx2);
            }
            if j + 1 < ny {
                rows.push(node);
                cols.push(node + 1);
                vals.push(-inv_dx2);
            }
            if i > 0 {
                rows.push(node);
                cols.push(node - ny);
                vals.push(-inv_dx2);
            }
            if i + 1 < nx {
                rows.push(node);
                cols.push(node + ny);
                vals.push(-inv_dx2);
            }
        }
    }
    CsrMatrix::from_triplets(n, n, &rows, &cols, &vals)
}
/// Conjugate Gradient solver for symmetric positive-definite systems `A x = b`.
///
/// Returns `(x, n_iter, residual_norm)`.
pub fn cg_solve(
    a: &CsrMatrix,
    b: &[f64],
    x0: &[f64],
    tol: f64,
    max_iter: usize,
) -> (Vec<f64>, u32, f64) {
    let n = b.len();
    assert_eq!(a.nrows, n);
    assert_eq!(a.ncols, n);
    let mut x = x0.to_vec();
    let ax0 = a.matvec(&x);
    let mut r: Vec<f64> = b.iter().zip(ax0.iter()).map(|(bi, ai)| bi - ai).collect();
    let mut p = r.clone();
    let mut rr: f64 = r.iter().map(|ri| ri * ri).sum();
    for iter in 0..max_iter {
        let res_norm = rr.sqrt();
        if res_norm < tol {
            return (x, iter as u32, res_norm);
        }
        let ap = a.matvec(&p);
        let pap: f64 = p.iter().zip(ap.iter()).map(|(pi, api)| pi * api).sum();
        if pap.abs() < f64::EPSILON {
            return (x, iter as u32, res_norm);
        }
        let alpha = rr / pap;
        for ((xi, ri), (&pi, &api)) in x.iter_mut().zip(r.iter_mut()).zip(p.iter().zip(ap.iter())) {
            *xi += alpha * pi;
            *ri -= alpha * api;
        }
        let rr_new: f64 = r.iter().map(|ri| ri * ri).sum();
        let beta = rr_new / rr;
        for (pi, &ri) in p.iter_mut().zip(r.iter()) {
            *pi = ri + beta * *pi;
        }
        rr = rr_new;
    }
    let res_norm = rr.sqrt();
    (x, max_iter as u32, res_norm)
}
/// Preconditioned Conjugate Gradient solver.
///
/// `precond[i]` = `1 / a_ii` (diagonal preconditioner).
/// Returns `(x, n_iter, residual_norm)`.
pub fn pcg_solve(
    a: &CsrMatrix,
    b: &[f64],
    precond: &[f64],
    tol: f64,
    max_iter: usize,
) -> (Vec<f64>, u32, f64) {
    let n = b.len();
    assert_eq!(a.nrows, n);
    assert_eq!(a.ncols, n);
    let mut x = vec![0.0f64; n];
    let ax0 = a.matvec(&x);
    let mut r: Vec<f64> = b.iter().zip(ax0.iter()).map(|(bi, ai)| bi - ai).collect();
    let mut z: Vec<f64> = r
        .iter()
        .zip(precond.iter())
        .map(|(ri, mi)| ri * mi)
        .collect();
    let mut p = z.clone();
    let mut rz: f64 = r.iter().zip(z.iter()).map(|(ri, zi)| ri * zi).sum();
    for iter in 0..max_iter {
        let res_norm: f64 = r.iter().map(|ri| ri * ri).sum::<f64>().sqrt();
        if res_norm < tol {
            return (x, iter as u32, res_norm);
        }
        let ap = a.matvec(&p);
        let pap: f64 = p.iter().zip(ap.iter()).map(|(pi, api)| pi * api).sum();
        if pap.abs() < f64::EPSILON {
            return (x, iter as u32, res_norm);
        }
        let alpha = rz / pap;
        for ((xi, ri), (&pi, &api)) in x.iter_mut().zip(r.iter_mut()).zip(p.iter().zip(ap.iter())) {
            *xi += alpha * pi;
            *ri -= alpha * api;
        }
        z = r
            .iter()
            .zip(precond.iter())
            .map(|(ri, mi)| ri * mi)
            .collect();
        let rz_new: f64 = r.iter().zip(z.iter()).map(|(ri, zi)| ri * zi).sum();
        let beta = rz_new / rz;
        for (pi, &zi) in p.iter_mut().zip(z.iter()) {
            *pi = zi + beta * *pi;
        }
        rz = rz_new;
    }
    let res_norm: f64 = r.iter().map(|ri| ri * ri).sum::<f64>().sqrt();
    (x, max_iter as u32, res_norm)
}
/// Successive Over-Relaxation (SOR) solver. Use `omega = 1` for Gauss-Seidel.
///
/// Returns `(x, n_iter)`.
pub fn sor_solve(
    a: &CsrMatrix,
    b: &[f64],
    omega: f64,
    tol: f64,
    max_iter: usize,
) -> (Vec<f64>, u32) {
    let n = b.len();
    assert_eq!(a.nrows, n);
    assert_eq!(a.ncols, n);
    let mut x = vec![0.0f64; n];
    for iter in 0..max_iter {
        let mut max_change = 0.0f64;
        for i in 0..n {
            let aii = a.get(i, i);
            if aii.abs() < f64::EPSILON {
                continue;
            }
            let mut sigma = 0.0;
            for k in a.row_ptr[i]..a.row_ptr[i + 1] {
                let j = a.col_idx[k];
                if j != i {
                    sigma += a.values[k] * x[j];
                }
            }
            let x_new = (b[i] - sigma) / aii;
            let delta = omega * (x_new - x[i]);
            max_change = max_change.max(delta.abs());
            x[i] += delta;
        }
        if max_change < tol {
            return (x, iter as u32 + 1);
        }
    }
    (x, max_iter as u32)
}
/// Jacobi iterative solver.
///
/// Returns `(x, n_iter)`.
pub fn jacobi_solve(a: &CsrMatrix, b: &[f64], tol: f64, max_iter: usize) -> (Vec<f64>, u32) {
    let n = b.len();
    assert_eq!(a.nrows, n);
    assert_eq!(a.ncols, n);
    let mut x = vec![0.0f64; n];
    for iter in 0..max_iter {
        let mut x_new = vec![0.0f64; n];
        let mut max_change = 0.0f64;
        for i in 0..n {
            let aii = a.get(i, i);
            if aii.abs() < f64::EPSILON {
                x_new[i] = x[i];
                continue;
            }
            let mut sigma = 0.0;
            for k in a.row_ptr[i]..a.row_ptr[i + 1] {
                let j = a.col_idx[k];
                if j != i {
                    sigma += a.values[k] * x[j];
                }
            }
            x_new[i] = (b[i] - sigma) / aii;
            max_change = max_change.max((x_new[i] - x[i]).abs());
        }
        x = x_new;
        if max_change < tol {
            return (x, iter as u32 + 1);
        }
    }
    (x, max_iter as u32)
}
/// Build a diagonal (Jacobi) preconditioner: returns `1 / a_ii` for each `i`.
///
/// Returns 1.0 where the diagonal is zero (to avoid division by zero).
pub fn diagonal_preconditioner(a: &CsrMatrix) -> Vec<f64> {
    let n = a.nrows.min(a.ncols);
    let mut m = vec![1.0f64; a.nrows];
    for (i, mi) in m.iter_mut().enumerate().take(n) {
        let d = a.get(i, i);
        if d.abs() > f64::EPSILON {
            *mi = 1.0 / d;
        }
    }
    m
}
/// Compute the dominant eigenvalue and eigenvector of a sparse matrix
/// using the power method.
///
/// Returns `(eigenvalue, eigenvector)` or `None` if not converged.
pub fn sparse_eig_power(a: &CsrMatrix, tol: f64, max_iter: usize) -> Option<(f64, Vec<f64>)> {
    let n = a.nrows;
    assert_eq!(n, a.ncols, "Matrix must be square");
    if n == 0 {
        return None;
    }
    let mut v: Vec<f64> = vec![1.0 / (n as f64).sqrt(); n];
    let mut eigenvalue = 0.0f64;
    for _ in 0..max_iter {
        let w = a.matvec(&v);
        let norm: f64 = w.iter().map(|x| x * x).sum::<f64>().sqrt();
        if norm < 1e-30 {
            return None;
        }
        let lambda_new: f64 = v.iter().zip(w.iter()).map(|(vi, wi)| vi * wi).sum();
        if (lambda_new - eigenvalue).abs() < tol {
            return Some((lambda_new, v));
        }
        eigenvalue = lambda_new;
        v = w.iter().map(|x| x / norm).collect();
    }
    Some((eigenvalue, v))
}
/// Compute the smallest eigenvalue of a sparse SPD matrix using inverse
/// iteration (shift-and-invert with the ILU preconditioner).
///
/// Returns `(eigenvalue, eigenvector)` or `None` if not converged.
pub fn sparse_eig_smallest(a: &CsrMatrix, tol: f64, max_iter: usize) -> Option<(f64, Vec<f64>)> {
    let n = a.nrows;
    assert_eq!(n, a.ncols, "Matrix must be square");
    if n == 0 {
        return None;
    }
    let lu = SparseLu::ilu0(a);
    let mut v: Vec<f64> = vec![1.0 / (n as f64).sqrt(); n];
    let mut eigenvalue = 0.0f64;
    for _ in 0..max_iter {
        let w = lu.solve(&v);
        let norm: f64 = w.iter().map(|x| x * x).sum::<f64>().sqrt();
        if norm < 1e-30 {
            return None;
        }
        let rq_inv: f64 = v.iter().zip(w.iter()).map(|(vi, wi)| vi * wi).sum();
        let lambda_new = if rq_inv.abs() < 1e-30 {
            0.0
        } else {
            1.0 / rq_inv
        };
        if (lambda_new - eigenvalue).abs() < tol {
            return Some((lambda_new, v));
        }
        eigenvalue = lambda_new;
        v = w.iter().map(|x| x / norm).collect();
    }
    Some((eigenvalue, v))
}
/// Estimate the spectral radius (largest singular value approximation) via
/// the power method applied to A^T A.
///
/// Returns the square root of the dominant eigenvalue of A^T A.
pub fn spectral_radius_estimate(a: &CsrMatrix, tol: f64, max_iter: usize) -> f64 {
    let at = a.transpose();
    let ata_matvec = |v: &[f64]| -> Vec<f64> {
        let w = a.matvec(v);
        at.matvec(&w)
    };
    let n = a.ncols;
    let mut v: Vec<f64> = vec![1.0 / (n as f64).sqrt(); n];
    let mut eigenvalue = 0.0f64;
    for _ in 0..max_iter {
        let w = ata_matvec(&v);
        let norm: f64 = w.iter().map(|x| x * x).sum::<f64>().sqrt();
        if norm < 1e-30 {
            break;
        }
        let lambda_new: f64 = v.iter().zip(w.iter()).map(|(vi, wi)| vi * wi).sum();
        if (lambda_new - eigenvalue).abs() < tol {
            return lambda_new.sqrt();
        }
        eigenvalue = lambda_new;
        v = w.iter().map(|x| x / norm).collect();
    }
    eigenvalue.sqrt()
}
/// Sparse-sparse matrix multiply: C = A * B.
///
/// A is `m × k`, B is `k × n`, result C is `m × n`.
/// Uses a dense accumulator row for each output row.
pub fn spgemm(a: &CsrMatrix, b: &CsrMatrix) -> CsrMatrix {
    assert_eq!(a.ncols, b.nrows, "spgemm: inner dimension mismatch");
    let m = a.nrows;
    let n = b.ncols;
    let mut c_row_ptr = vec![0usize; m + 1];
    let mut c_col_idx: Vec<usize> = Vec::new();
    let mut c_values: Vec<f64> = Vec::new();
    let mut accum = vec![0.0f64; n];
    let mut marker = vec![usize::MAX; n];
    for i in 0..m {
        let mut touched: Vec<usize> = Vec::new();
        for ka in a.row_ptr[i]..a.row_ptr[i + 1] {
            let p = a.col_idx[ka];
            let a_ip = a.values[ka];
            for kb in b.row_ptr[p]..b.row_ptr[p + 1] {
                let j = b.col_idx[kb];
                accum[j] += a_ip * b.values[kb];
                if marker[j] != i {
                    marker[j] = i;
                    touched.push(j);
                }
            }
        }
        touched.sort_unstable();
        for j in &touched {
            if accum[*j].abs() > 0.0 {
                c_col_idx.push(*j);
                c_values.push(accum[*j]);
            }
            accum[*j] = 0.0;
        }
        c_row_ptr[i + 1] = c_col_idx.len();
    }
    CsrMatrix {
        nrows: m,
        ncols: n,
        row_ptr: c_row_ptr,
        col_idx: c_col_idx,
        values: c_values,
    }
}
/// Convert COO (coordinate format) triplets directly to CSR.
///
/// Handles unsorted input; duplicate (row, col) entries are summed.
/// This is a standalone function complementing `CsrMatrix::from_triplets`.
pub fn coo_to_csr(
    nrows: usize,
    ncols: usize,
    rows: &[usize],
    cols: &[usize],
    vals: &[f64],
) -> CsrMatrix {
    CsrMatrix::from_triplets(nrows, ncols, rows, cols, vals)
}
/// Arnoldi iteration: build an orthonormal Krylov basis for `A`.
///
/// Starting from initial vector `b0`, computes `k` basis vectors
/// `Q = [q_0, q_1, ..., q_{k-1}]` (each of length `n`) and the
/// upper Hessenberg matrix `H` (stored row-major as `Vec<Vec`f64`>`
/// of size `(k+1) × k`).
///
/// The Arnoldi relation holds: A Q_k = Q_{k+1} H_{k+1,k}.
///
/// # Arguments
/// * `a`  – input sparse matrix (n × n)
/// * `b0` – starting vector (length n)
/// * `k`  – number of Arnoldi steps to perform
///
/// # Returns
/// `(Q, H)` where `Q[i]` is the i-th Krylov basis vector and
/// `H[i][j]` is the Hessenberg coefficient.
pub fn arnoldi(a: &CsrMatrix, b0: &[f64], k: usize) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let n = a.nrows;
    assert_eq!(
        b0.len(),
        n,
        "arnoldi: b0 length must equal matrix dimension"
    );
    let norm0: f64 = b0.iter().map(|x| x * x).sum::<f64>().sqrt();
    assert!(norm0 > 1e-30, "arnoldi: starting vector is zero");
    let mut q: Vec<Vec<f64>> = Vec::with_capacity(k + 1);
    let mut h: Vec<Vec<f64>> = vec![vec![0.0; k]; k + 1];
    q.push(b0.iter().map(|x| x / norm0).collect());
    for j in 0..k {
        if j >= q.len() {
            break;
        }
        let mut w = a.matvec(&q[j]);
        for i in 0..=j {
            let hij: f64 = w.iter().zip(q[i].iter()).map(|(wi, qi)| wi * qi).sum();
            h[i][j] = hij;
            for (wl, qi) in w.iter_mut().zip(q[i].iter()) {
                *wl -= hij * qi;
            }
        }
        let beta: f64 = w.iter().map(|x| x * x).sum::<f64>().sqrt();
        h[j + 1][j] = beta;
        if beta > 1e-14 {
            q.push(w.iter().map(|x| x / beta).collect());
        } else {
            break;
        }
    }
    let q_len = q.len();
    h.truncate(q_len);
    (q, h)
}
/// Diagonal (Jacobi) scaling: return D^{-1} A where D = diag(A).
///
/// Rows with zero diagonal are left unchanged (scale factor = 1).
/// This is the simplest preconditioner and is used to equilibrate the
/// system before iterative solves.
pub fn jacobi_scale(a: &CsrMatrix) -> CsrMatrix {
    let diag = a.diagonal();
    let mut scaled = CsrMatrix {
        nrows: a.nrows,
        ncols: a.ncols,
        row_ptr: a.row_ptr.clone(),
        col_idx: a.col_idx.clone(),
        values: a.values.clone(),
    };
    for (i, &d) in diag.iter().enumerate().take(a.nrows) {
        let scale = if d.abs() > 1e-30 { 1.0 / d } else { 1.0 };
        for k in a.row_ptr[i]..a.row_ptr[i + 1] {
            scaled.values[k] *= scale;
        }
    }
    scaled
}
/// Biconjugate Gradient Stabilised (BICGStab) solver for non-symmetric systems A x = b.
///
/// Returns `(solution, converged, iterations)`.
pub fn bicgstab_solve(
    a: &CsrMatrix,
    b: &[f64],
    x0: Option<&[f64]>,
    tol: f64,
    max_iter: usize,
) -> (Vec<f64>, bool, usize) {
    let n = a.nrows;
    assert_eq!(b.len(), n);
    let mut x: Vec<f64> = match x0 {
        Some(v) => v.to_vec(),
        None => vec![0.0; n],
    };
    let ax = a.matvec(&x);
    let mut r: Vec<f64> = b.iter().zip(ax.iter()).map(|(bi, ai)| bi - ai).collect();
    let r_hat = r.clone();
    let mut rho_prev = 1.0f64;
    let mut alpha = 1.0f64;
    let mut omega = 1.0f64;
    let mut p = vec![0.0f64; n];
    let mut v = vec![0.0f64; n];
    let b_norm: f64 = b.iter().map(|x| x * x).sum::<f64>().sqrt().max(1e-30);
    for iter in 0..max_iter {
        let rho = r_hat.iter().zip(r.iter()).map(|(a, b)| a * b).sum::<f64>();
        if rho.abs() < 1e-300 {
            break;
        }
        if iter == 0 {
            p = r.clone();
        } else {
            let beta = (rho / rho_prev) * (alpha / omega);
            for ((pi, &ri), &vi) in p.iter_mut().zip(r.iter()).zip(v.iter()) {
                *pi = ri + beta * (*pi - omega * vi);
            }
        }
        v = a.matvec(&p);
        let r_hat_v: f64 = r_hat.iter().zip(v.iter()).map(|(a, b)| a * b).sum();
        if r_hat_v.abs() < 1e-300 {
            break;
        }
        alpha = rho / r_hat_v;
        let s: Vec<f64> = r
            .iter()
            .zip(v.iter())
            .map(|(ri, vi)| ri - alpha * vi)
            .collect();
        let s_norm: f64 = s.iter().map(|x| x * x).sum::<f64>().sqrt();
        if s_norm / b_norm < tol {
            for (xi, &pi) in x.iter_mut().zip(p.iter()) {
                *xi += alpha * pi;
            }
            return (x, true, iter + 1);
        }
        let t = a.matvec(&s);
        let t_dot_s: f64 = t.iter().zip(s.iter()).map(|(a, b)| a * b).sum();
        let t_dot_t: f64 = t.iter().map(|x| x * x).sum();
        omega = if t_dot_t < 1e-300 {
            0.0
        } else {
            t_dot_s / t_dot_t
        };
        for ((xi, &pi), &si) in x.iter_mut().zip(p.iter()).zip(s.iter()) {
            *xi += alpha * pi + omega * si;
        }
        for (ri, (&si, &ti)) in r.iter_mut().zip(s.iter().zip(t.iter())) {
            *ri = si - omega * ti;
        }
        let r_norm: f64 = r.iter().map(|x| x * x).sum::<f64>().sqrt();
        if r_norm / b_norm < tol {
            return (x, true, iter + 1);
        }
        if omega.abs() < 1e-300 {
            break;
        }
        rho_prev = rho;
    }
    let r_final = {
        let ax2 = a.matvec(&x);
        let res: Vec<f64> = b.iter().zip(ax2.iter()).map(|(bi, ai)| bi - ai).collect();
        res.iter().map(|x| x * x).sum::<f64>().sqrt()
    };
    (x, r_final / b_norm < tol, max_iter)
}
/// Minimum Residual (MINRES) method for symmetric (possibly indefinite) systems A x = b.
///
/// Returns `(solution, converged, iterations)`.
pub fn minres_solve(
    a: &CsrMatrix,
    b: &[f64],
    tol: f64,
    max_iter: usize,
) -> (Vec<f64>, bool, usize) {
    let n = a.nrows;
    assert_eq!(b.len(), n);
    let mut x = vec![0.0f64; n];
    let mut r = b.to_vec();
    let b_norm: f64 = b.iter().map(|x| x * x).sum::<f64>().sqrt().max(1e-30);
    let mut v_old = vec![0.0f64; n];
    let mut beta: f64 = r.iter().map(|x| x * x).sum::<f64>().sqrt();
    if beta < 1e-30 {
        return (x, true, 0);
    }
    let mut v = r.iter().map(|x| x / beta).collect::<Vec<f64>>();
    let mut d = vec![0.0f64; n];
    let mut d_old = vec![0.0f64; n];
    let mut phi_bar = beta;
    let mut c_old = 1.0f64;
    let mut s_old = 0.0f64;
    let mut c = 1.0f64;
    let mut s = 0.0f64;
    for iter in 0..max_iter {
        let av = a.matvec(&v);
        let alpha: f64 = v.iter().zip(av.iter()).map(|(vi, avi)| vi * avi).sum();
        let mut v_new: Vec<f64> = av
            .iter()
            .enumerate()
            .map(|(i, &avi)| avi - alpha * v[i] - beta * v_old[i])
            .collect();
        let beta_new: f64 = v_new.iter().map(|x| x * x).sum::<f64>().sqrt();
        if beta_new > 1e-30 {
            for vi in &mut v_new {
                *vi /= beta_new;
            }
        }
        let eps = s_old * beta;
        let delta = c_old * c * beta + s * alpha;
        let gamma_bar = -s_old * s * beta + c * alpha - c_old * delta * s_old;
        let rhs = phi_bar * c;
        let phi = phi_bar * s;
        let _ = eps;
        let _ = delta;
        let _ = gamma_bar;
        let _ = rhs;
        let _ = phi;
        let gamma: f64 = (c * alpha + s * beta_new).hypot(beta_new);
        if gamma < 1e-30 {
            break;
        }
        let c_new = (c * alpha + s * beta_new) / gamma;
        let s_new = beta_new / gamma;
        let d_new: Vec<f64> = v
            .iter()
            .enumerate()
            .map(|(i, &vi)| (vi - c * d[i] - s * d_old[i]) / gamma)
            .collect();
        let eta = c_new * phi_bar;
        for (xi, &dni) in x.iter_mut().zip(d_new.iter()) {
            *xi += eta * dni;
        }
        phi_bar *= s_new;
        if phi_bar / b_norm < tol {
            return (x, true, iter + 1);
        }
        let _ = c_old;
        let _ = s_old;
        c_old = c;
        s_old = s;
        c = c_new;
        s = s_new;
        d_old = d;
        d = d_new;
        v_old = v;
        v = v_new;
        beta = beta_new;
        r = b
            .iter()
            .zip(a.matvec(&x).iter())
            .map(|(bi, ai)| bi - ai)
            .collect();
        let r_norm: f64 = r.iter().map(|x| x * x).sum::<f64>().sqrt();
        if r_norm / b_norm < tol {
            return (x, true, iter + 1);
        }
    }
    (x, false, max_iter)
}
/// Gauss-Seidel iterative solver for A x = b.
///
/// Requires A to have non-zero diagonal entries.
/// Returns `(solution, iterations_performed)`.
pub fn gauss_seidel_solve(
    a: &CsrMatrix,
    b: &[f64],
    tol: f64,
    max_iter: usize,
) -> (Vec<f64>, usize) {
    let n = a.nrows;
    assert_eq!(b.len(), n);
    let mut x = vec![0.0f64; n];
    for iter in 0..max_iter {
        let mut max_change = 0.0f64;
        for i in 0..n {
            let diag = a.get(i, i);
            assert!(diag.abs() > 1e-30, "Gauss-Seidel: zero diagonal at row {i}");
            let row_sum: f64 = (a.row_ptr[i]..a.row_ptr[i + 1])
                .map(|k| {
                    let j = a.col_idx[k];
                    if j != i { a.values[k] * x[j] } else { 0.0 }
                })
                .sum();
            let x_new = (b[i] - row_sum) / diag;
            max_change = max_change.max((x_new - x[i]).abs());
            x[i] = x_new;
        }
        if max_change < tol {
            return (x, iter + 1);
        }
    }
    (x, max_iter)
}
/// Forward substitution: solve L x = b where L is lower triangular (CSR).
pub fn forward_substitution(l: &CsrMatrix, b: &[f64]) -> Vec<f64> {
    let n = l.nrows;
    assert_eq!(b.len(), n);
    let mut x = vec![0.0f64; n];
    for i in 0..n {
        let mut sum = b[i];
        let diag = l.get(i, i);
        for k in l.row_ptr[i]..l.row_ptr[i + 1] {
            let j = l.col_idx[k];
            if j < i {
                sum -= l.values[k] * x[j];
            }
        }
        assert!(
            diag.abs() > 1e-30,
            "forward_substitution: zero diagonal at row {i}"
        );
        x[i] = sum / diag;
    }
    x
}
/// Backward substitution: solve U x = b where U is upper triangular (CSR).
pub fn backward_substitution(u: &CsrMatrix, b: &[f64]) -> Vec<f64> {
    let n = u.nrows;
    assert_eq!(b.len(), n);
    let mut x = vec![0.0f64; n];
    for i in (0..n).rev() {
        let mut sum = b[i];
        let diag = u.get(i, i);
        for k in u.row_ptr[i]..u.row_ptr[i + 1] {
            let j = u.col_idx[k];
            if j > i {
                sum -= u.values[k] * x[j];
            }
        }
        assert!(
            diag.abs() > 1e-30,
            "backward_substitution: zero diagonal at row {i}"
        );
        x[i] = sum / diag;
    }
    x
}
/// Build a CSR matrix from a banded stencil.
///
/// `diagonals[k]` is the value at offset `offsets[k]` (0 = main diagonal,
/// positive = superdiagonal, negative = subdiagonal).
pub fn banded_matrix(n: usize, offsets: &[i64], diagonals: &[f64]) -> CsrMatrix {
    assert_eq!(offsets.len(), diagonals.len());
    let mut rows = Vec::new();
    let mut cols = Vec::new();
    let mut vals = Vec::new();
    for (k, &off) in offsets.iter().enumerate() {
        for i in 0..n {
            let j = i as i64 + off;
            if j >= 0 && (j as usize) < n {
                rows.push(i);
                cols.push(j as usize);
                vals.push(diagonals[k]);
            }
        }
    }
    CsrMatrix::from_triplets(n, n, &rows, &cols, &vals)
}
/// Row-scale a CSR matrix: D_r A where D_r = diag(scales).
pub fn row_scale(a: &CsrMatrix, scales: &[f64]) -> CsrMatrix {
    assert_eq!(scales.len(), a.nrows);
    let values: Vec<f64> = (0..a.nrows)
        .flat_map(|i| {
            let s = scales[i];
            (a.row_ptr[i]..a.row_ptr[i + 1]).map(move |k| a.values[k] * s)
        })
        .collect();
    CsrMatrix {
        nrows: a.nrows,
        ncols: a.ncols,
        row_ptr: a.row_ptr.clone(),
        col_idx: a.col_idx.clone(),
        values,
    }
}
/// Column-scale a CSR matrix: A D_c where D_c = diag(scales).
pub fn col_scale(a: &CsrMatrix, scales: &[f64]) -> CsrMatrix {
    assert_eq!(scales.len(), a.ncols);
    let values: Vec<f64> = a
        .col_idx
        .iter()
        .zip(a.values.iter())
        .map(|(&j, &v)| v * scales[j])
        .collect();
    CsrMatrix {
        nrows: a.nrows,
        ncols: a.ncols,
        row_ptr: a.row_ptr.clone(),
        col_idx: a.col_idx.clone(),
        values,
    }
}
/// Two-sided equilibration scaling: compute row/col scale vectors so that
/// the scaled matrix has unit diagonal (i.e., sqrt(|A_ii|) scaling).
pub fn equilibration_scales(a: &CsrMatrix) -> (Vec<f64>, Vec<f64>) {
    let n = a.nrows.max(a.ncols);
    let mut row_scales = vec![1.0f64; a.nrows];
    let mut col_scales = vec![1.0f64; a.ncols];
    for (i, rs) in row_scales.iter_mut().enumerate() {
        let d = a.get(i, i).abs();
        if d > 1e-30 {
            *rs = 1.0 / d.sqrt();
        }
    }
    for (j, cs) in col_scales.iter_mut().enumerate() {
        let d = a.get(j.min(a.nrows - 1), j).abs();
        if d > 1e-30 {
            *cs = 1.0 / d.sqrt();
        }
    }
    let _ = n;
    (row_scales, col_scales)
}
/// Solve A X = B for multiple right-hand sides simultaneously using CG.
///
/// Each column of `b_cols` is a separate RHS; returns a `Vec<Vec`f64`>` of solutions.
pub fn multi_rhs_cg(
    a: &CsrMatrix,
    b_cols: &[Vec<f64>],
    tol: f64,
    max_iter: usize,
) -> Vec<Vec<f64>> {
    b_cols
        .iter()
        .map(|b| {
            let x0 = vec![0.0f64; b.len()];
            let (x, _iters, _res) = cg_solve(a, b, &x0, tol, max_iter);
            x
        })
        .collect()
}
/// Solve the saddle-point system:
///   \[ A   B^T \] \[ x \]   \[ f \]
///   \[ B   0   \] \[ p \] = \[ g \]
///
/// Uses a Schur complement approach with CG on the pressure sub-system.
///
/// # Arguments
/// * `a`        – (n×n) SPD matrix
/// * `bt`       – (n×m) matrix B^T
/// * `b`        – (m×n) matrix B
/// * `f`        – RHS vector of length n
/// * `g`        – RHS vector of length m
/// * `tol`      – solver tolerance
/// * `max_iter` – max iterations
pub fn saddle_point_solve(
    a: &CsrMatrix,
    bt: &CsrMatrix,
    b: &CsrMatrix,
    f: &[f64],
    g: &[f64],
    tol: f64,
    max_iter: usize,
) -> (Vec<f64>, Vec<f64>) {
    let n = a.nrows;
    let m = b.nrows;
    assert_eq!(f.len(), n);
    assert_eq!(g.len(), m);
    let x0_init = vec![0.0f64; n];
    let (x0, _iters0, _res0) = cg_solve(a, f, &x0_init, tol, max_iter);
    let b_x0 = b.matvec(&x0);
    let schur_rhs: Vec<f64> = g
        .iter()
        .zip(b_x0.iter())
        .map(|(gi, bxi)| gi - bxi)
        .collect();
    let s_diag: Vec<f64> = (0..m)
        .map(|i| {
            let mut e = vec![0.0f64; m];
            e[i] = 1.0;
            let bt_e = bt.matvec(&e);
            let x0_schur = vec![0.0f64; n];
            let (a_inv_bt_e, _, _) = cg_solve(a, &bt_e, &x0_schur, tol * 0.1, max_iter);
            let b_a_inv = b.matvec(&a_inv_bt_e);
            b_a_inv[i]
        })
        .collect();
    let p: Vec<f64> = schur_rhs
        .iter()
        .zip(s_diag.iter())
        .map(|(r, s)| if s.abs() > 1e-30 { r / s } else { 0.0 })
        .collect();
    let bt_p = bt.matvec(&p);
    let f_corr: Vec<f64> = f.iter().zip(bt_p.iter()).map(|(fi, bi)| fi - bi).collect();
    let x0_final = vec![0.0f64; n];
    let (x, _, _) = cg_solve(a, &f_corr, &x0_final, tol, max_iter);
    (x, p)
}
/// Lanczos iteration: build a tridiagonal approximation of a symmetric matrix A.
///
/// Returns `(alphas, betas, Q)` where:
/// - `alphas[k]` is the k-th diagonal of the tridiagonal matrix,
/// - `betas[k]` is the k-th subdiagonal (length k steps),
/// - `Q[k]` is the k-th Lanczos vector (orthonormal).
pub fn lanczos(a: &CsrMatrix, b0: &[f64], k: usize) -> (Vec<f64>, Vec<f64>, Vec<Vec<f64>>) {
    let n = a.nrows;
    assert_eq!(b0.len(), n);
    let norm0: f64 = b0.iter().map(|x| x * x).sum::<f64>().sqrt();
    assert!(norm0 > 1e-30, "lanczos: starting vector is zero");
    let mut q: Vec<Vec<f64>> = Vec::with_capacity(k + 1);
    q.push(b0.iter().map(|x| x / norm0).collect());
    let mut alphas = Vec::with_capacity(k);
    let mut betas = Vec::with_capacity(k);
    let mut q_prev = vec![0.0f64; n];
    for j in 0..k {
        let mut w = a.matvec(&q[j]);
        let alpha: f64 = q[j].iter().zip(w.iter()).map(|(qi, wi)| qi * wi).sum();
        alphas.push(alpha);
        if j > 0 {
            let beta_prev = betas[j - 1];
            for ((wi, &qji), &qpi) in w.iter_mut().zip(q[j].iter()).zip(q_prev.iter()) {
                *wi -= alpha * qji + beta_prev * qpi;
            }
        } else {
            for (wi, &qji) in w.iter_mut().zip(q[j].iter()) {
                *wi -= alpha * qji;
            }
        }
        let beta: f64 = w.iter().map(|x| x * x).sum::<f64>().sqrt();
        betas.push(beta);
        if beta < 1e-30 || j + 1 >= k {
            break;
        }
        q_prev = q[j].clone();
        q.push(w.iter().map(|x| x / beta).collect());
    }
    (alphas, betas, q)
}
/// Compute Ritz values from a Lanczos tridiagonal matrix via QR iteration.
///
/// Given `alphas` (diagonal) and `betas` (subdiagonal), returns approximate eigenvalues.
pub fn ritz_values(alphas: &[f64], betas: &[f64]) -> Vec<f64> {
    let m = alphas.len();
    if m == 0 {
        return Vec::new();
    }
    let mut diag = alphas.to_vec();
    let mut off = betas[..betas.len().min(m - 1)].to_vec();
    for _ in 0..500 {
        let shift = diag[m - 1];
        for d in &mut diag {
            *d -= shift;
        }
        let mut c = vec![1.0f64; m];
        let mut s_vec = vec![0.0f64; m];
        let mut d_new = diag.clone();
        let mut off_new = off.clone();
        for i in 0..m - 1 {
            let r = (diag[i] * diag[i] + off[i] * off[i]).sqrt();
            if r < 1e-30 {
                continue;
            }
            c[i] = diag[i] / r;
            s_vec[i] = off[i] / r;
            let new_d = c[i] * diag[i] + s_vec[i] * off[i];
            d_new[i] = new_d;
            if i + 1 < m {
                off_new[i] = c[i] * off[i] - s_vec[i] * diag[i];
                d_new[i + 1] = -s_vec[i] * off[i] + c[i] * diag[i + 1];
            }
        }
        diag = d_new;
        off = off_new;
        for d in &mut diag {
            *d += shift;
        }
        let off_norm: f64 = off.iter().map(|x| x * x).sum::<f64>().sqrt();
        if off_norm < 1e-10 {
            break;
        }
    }
    diag
}
/// Invert a small dense matrix (size × size) stored row-major using Gauss-Jordan.
/// Returns `None` if singular.
pub(super) fn invert_small_dense(m: &[f64], size: usize) -> Option<Vec<f64>> {
    let mut a = m.to_vec();
    let mut inv = vec![0.0f64; size * size];
    for i in 0..size {
        inv[i * size + i] = 1.0;
    }
    for col in 0..size {
        let mut pivot_row = col;
        let mut max_val = a[col * size + col].abs();
        for row in (col + 1)..size {
            if a[row * size + col].abs() > max_val {
                max_val = a[row * size + col].abs();
                pivot_row = row;
            }
        }
        if max_val < 1e-14 {
            return None;
        }
        for j in 0..size {
            a.swap(col * size + j, pivot_row * size + j);
            inv.swap(col * size + j, pivot_row * size + j);
        }
        let piv = a[col * size + col];
        for j in 0..size {
            a[col * size + j] /= piv;
            inv[col * size + j] /= piv;
        }
        for row in 0..size {
            if row == col {
                continue;
            }
            let factor = a[row * size + col];
            for j in 0..size {
                let av = a[col * size + j];
                a[row * size + j] -= factor * av;
                let iv = inv[col * size + j];
                inv[row * size + j] -= factor * iv;
            }
        }
    }
    Some(inv)
}
/// Extract the lower triangular part of a sparse matrix (including diagonal).
pub fn lower_triangular(a: &CsrMatrix) -> CsrMatrix {
    let mut rows = Vec::new();
    let mut cols = Vec::new();
    let mut vals = Vec::new();
    for i in 0..a.nrows {
        for k in a.row_ptr[i]..a.row_ptr[i + 1] {
            let j = a.col_idx[k];
            if j <= i {
                rows.push(i);
                cols.push(j);
                vals.push(a.values[k]);
            }
        }
    }
    CsrMatrix::from_triplets(a.nrows, a.ncols, &rows, &cols, &vals)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::sparse::CscMatrix;
    use crate::sparse::IncompleteCholesky;
    use crate::sparse::SparseQr;
    use crate::sparse::SsorPreconditioner;
    pub(super) const TOL: f64 = 1e-10;
    #[test]
    fn test_new_empty() {
        let m = CsrMatrix::new(3, 4);
        assert_eq!(m.nrows, 3);
        assert_eq!(m.ncols, 4);
        assert_eq!(m.nnz(), 0);
        assert_eq!(m.row_ptr.len(), 4);
    }
    #[test]
    fn test_from_triplets_basic() {
        let rows = [0, 1, 1, 2];
        let cols = [0, 0, 2, 2];
        let vals = [1.0, 2.0, 3.0, 4.0];
        let m = CsrMatrix::from_triplets(3, 3, &rows, &cols, &vals);
        assert_eq!(m.nnz(), 4);
        assert!((m.get(0, 0) - 1.0).abs() < TOL);
        assert!((m.get(1, 0) - 2.0).abs() < TOL);
        assert!((m.get(1, 2) - 3.0).abs() < TOL);
        assert!((m.get(2, 2) - 4.0).abs() < TOL);
        assert!((m.get(0, 1)).abs() < TOL);
    }
    #[test]
    fn test_from_triplets_duplicate_sum() {
        let rows = [0, 0];
        let cols = [0, 0];
        let vals = [1.5, 2.5];
        let m = CsrMatrix::from_triplets(2, 2, &rows, &cols, &vals);
        assert_eq!(m.nnz(), 1);
        assert!((m.get(0, 0) - 4.0).abs() < TOL);
    }
    #[test]
    fn test_get_zero_for_missing() {
        let m = CsrMatrix::new(4, 4);
        assert_eq!(m.get(2, 3), 0.0);
    }
    #[test]
    fn test_matvec_identity() {
        let id = identity_csr(4);
        let x = vec![1.0, 2.0, 3.0, 4.0];
        let y = id.matvec(&x);
        for (yi, xi) in y.iter().zip(x.iter()) {
            assert!((yi - xi).abs() < TOL);
        }
    }
    #[test]
    fn test_matvec_simple() {
        let m = CsrMatrix::from_triplets(2, 2, &[0, 0, 1], &[0, 1, 1], &[2.0, 1.0, 3.0]);
        let y = m.matvec(&[1.0, 2.0]);
        assert!((y[0] - 4.0).abs() < TOL);
        assert!((y[1] - 6.0).abs() < TOL);
    }
    #[test]
    fn test_add_matrices() {
        let a = identity_csr(3);
        let b = identity_csr(3);
        let c = a.add(&b);
        for i in 0..3 {
            assert!((c.get(i, i) - 2.0).abs() < TOL);
        }
    }
    #[test]
    fn test_add_overlap() {
        let r = [0usize, 0];
        let c = [0usize, 0];
        let v1 = [1.0, 0.0];
        let v2 = [0.0, 1.0];
        let a = CsrMatrix::from_triplets(2, 2, &r[..1], &c[..1], &v1[..1]);
        let b = CsrMatrix::from_triplets(2, 2, &r[..1], &c[..1], &v2[..1]);
        let s = a.add(&b);
        assert!((s.get(0, 0) - 1.0).abs() < TOL);
    }
    #[test]
    fn test_scale() {
        let id = identity_csr(3);
        let scaled = id.scale(5.0);
        for i in 0..3 {
            assert!((scaled.get(i, i) - 5.0).abs() < TOL);
        }
    }
    #[test]
    fn test_scale_zero() {
        let id = identity_csr(3);
        let z = id.scale(0.0);
        for i in 0..3 {
            assert!(z.get(i, i).abs() < TOL);
        }
    }
    #[test]
    fn test_transpose_identity() {
        let id = identity_csr(4);
        let t = id.transpose();
        let d = t.to_dense();
        let orig = id.to_dense();
        assert_eq!(d, orig);
    }
    #[test]
    fn test_transpose_rectangular() {
        let m = CsrMatrix::from_triplets(2, 3, &[0, 0, 1], &[0, 2, 1], &[1.0, 2.0, 3.0]);
        let t = m.transpose();
        assert_eq!(t.nrows, 3);
        assert_eq!(t.ncols, 2);
        assert!((t.get(0, 0) - 1.0).abs() < TOL);
        assert!((t.get(2, 0) - 2.0).abs() < TOL);
        assert!((t.get(1, 1) - 3.0).abs() < TOL);
    }
    #[test]
    fn test_diagonal_identity() {
        let id = identity_csr(5);
        let d = id.diagonal();
        assert_eq!(d.len(), 5);
        for v in d {
            assert!((v - 1.0).abs() < TOL);
        }
    }
    #[test]
    fn test_diagonal_general() {
        let m = CsrMatrix::from_triplets(3, 3, &[0, 1, 2, 0], &[0, 1, 2, 2], &[7.0, 8.0, 9.0, 1.0]);
        let d = m.diagonal();
        assert!((d[0] - 7.0).abs() < TOL);
        assert!((d[1] - 8.0).abs() < TOL);
        assert!((d[2] - 9.0).abs() < TOL);
    }
    #[test]
    fn test_to_dense_identity() {
        let id = identity_csr(3);
        let d = id.to_dense();
        for (i, row) in d.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((val - expected).abs() < TOL);
            }
        }
    }
    #[test]
    fn test_from_dense_round_trip() {
        let dense = vec![
            vec![1.0, 0.0, 2.0],
            vec![0.0, 3.0, 0.0],
            vec![4.0, 0.0, 5.0],
        ];
        let m = CsrMatrix::from_dense(&dense);
        let back = m.to_dense();
        for i in 0..3 {
            for j in 0..3 {
                assert!((back[i][j] - dense[i][j]).abs() < TOL);
            }
        }
    }
    #[test]
    fn test_is_symmetric_identity() {
        assert!(identity_csr(5).is_symmetric(TOL));
    }
    #[test]
    fn test_is_symmetric_laplacian() {
        let lap = laplacian_1d(10, 0.1);
        assert!(lap.is_symmetric(TOL));
    }
    #[test]
    fn test_is_not_symmetric() {
        let m = CsrMatrix::from_triplets(2, 2, &[0], &[1], &[1.0]);
        assert!(!m.is_symmetric(TOL));
    }
    #[test]
    fn test_identity_csr() {
        let id = identity_csr(6);
        assert_eq!(id.nnz(), 6);
        let y = id.matvec(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        for (i, yi) in y.iter().enumerate() {
            assert!((yi - (i + 1) as f64).abs() < TOL);
        }
    }
    #[test]
    fn test_laplacian_1d_size() {
        let n = 8;
        let lap = laplacian_1d(n, 0.1);
        assert_eq!(lap.nrows, n);
        assert_eq!(lap.ncols, n);
        assert!(lap.is_symmetric(TOL));
    }
    #[test]
    fn test_laplacian_1d_stencil() {
        let dx = 0.5;
        let lap = laplacian_1d(3, dx);
        let inv = 1.0 / (dx * dx);
        assert!((lap.get(0, 0) - 2.0 * inv).abs() < TOL);
        assert!((lap.get(0, 1) + inv).abs() < TOL);
        assert!((lap.get(1, 0) + inv).abs() < TOL);
        assert!((lap.get(1, 1) - 2.0 * inv).abs() < TOL);
    }
    #[test]
    fn test_laplacian_2d_size() {
        let (nx, ny) = (4, 5);
        let lap = laplacian_2d(nx, ny, 0.1);
        assert_eq!(lap.nrows, nx * ny);
        assert_eq!(lap.ncols, nx * ny);
        assert!(lap.is_symmetric(TOL));
    }
    #[test]
    fn test_laplacian_2d_diagonal() {
        let lap = laplacian_2d(3, 3, 1.0);
        let d = lap.diagonal();
        for v in d {
            assert!((v - 4.0).abs() < TOL);
        }
    }
    #[test]
    fn test_cg_solve_simple() {
        let a = CsrMatrix::from_triplets(2, 2, &[0, 0, 1, 1], &[0, 1, 0, 1], &[4.0, 1.0, 1.0, 3.0]);
        let b = [1.0, 2.0];
        let x0 = [0.0, 0.0];
        let (x, _, res) = cg_solve(&a, &b, &x0, 1e-10, 100);
        let ax = a.matvec(&x);
        assert!((ax[0] - b[0]).abs() < 1e-8);
        assert!((ax[1] - b[1]).abs() < 1e-8);
        assert!(res < 1e-8);
    }
    #[test]
    fn test_cg_solve_identity() {
        let n = 5;
        let id = identity_csr(n);
        let b: Vec<f64> = (1..=n).map(|i| i as f64).collect();
        let x0 = vec![0.0; n];
        let (x, _, _) = cg_solve(&id, &b, &x0, 1e-12, 100);
        for (xi, bi) in x.iter().zip(b.iter()) {
            assert!((xi - bi).abs() < 1e-10);
        }
    }
    #[test]
    fn test_cg_solve_laplacian() {
        let n = 10;
        let lap = laplacian_1d(n, 1.0 / (n + 1) as f64);
        let b = vec![1.0; n];
        let x0 = vec![0.0; n];
        let (x, _, res) = cg_solve(&lap, &b, &x0, 1e-10, 1000);
        let ax = lap.matvec(&x);
        let err: f64 = ax
            .iter()
            .zip(b.iter())
            .map(|(ai, bi)| (ai - bi).powi(2))
            .sum::<f64>()
            .sqrt();
        assert!(err < 1e-8, "residual={res}, err={err}");
    }
    #[test]
    fn test_jacobi_solve_diagonal() {
        let a = CsrMatrix::from_triplets(3, 3, &[0, 1, 2], &[0, 1, 2], &[2.0, 4.0, 8.0]);
        let b = [4.0, 8.0, 16.0];
        let (x, _) = jacobi_solve(&a, &b, 1e-12, 100);
        assert!((x[0] - 2.0).abs() < TOL);
        assert!((x[1] - 2.0).abs() < TOL);
        assert!((x[2] - 2.0).abs() < TOL);
    }
    #[test]
    fn test_jacobi_solve_convergence() {
        let a =
            CsrMatrix::from_triplets(2, 2, &[0, 0, 1, 1], &[0, 1, 0, 1], &[10.0, 1.0, 1.0, 10.0]);
        let b = [11.0, 11.0];
        let (x, _) = jacobi_solve(&a, &b, 1e-12, 1000);
        assert!((x[0] - 1.0).abs() < 1e-9);
        assert!((x[1] - 1.0).abs() < 1e-9);
    }
    #[test]
    fn test_sor_gauss_seidel() {
        let a = CsrMatrix::from_triplets(2, 2, &[0, 0, 1, 1], &[0, 1, 0, 1], &[4.0, 1.0, 1.0, 3.0]);
        let b = [1.0, 2.0];
        let (x, _) = sor_solve(&a, &b, 1.0, 1e-12, 10000);
        let ax = a.matvec(&x);
        assert!((ax[0] - b[0]).abs() < 1e-9);
        assert!((ax[1] - b[1]).abs() < 1e-9);
    }
    #[test]
    fn test_sor_omega() {
        let a = CsrMatrix::from_triplets(
            3,
            3,
            &[0, 1, 2, 0, 1, 1, 2],
            &[0, 1, 2, 1, 0, 2, 1],
            &[4.0, 4.0, 4.0, -1.0, -1.0, -1.0, -1.0],
        );
        let b = [3.0, 2.0, 3.0];
        let (x, _) = sor_solve(&a, &b, 1.2, 1e-10, 10000);
        let ax = a.matvec(&x);
        for i in 0..3 {
            assert!(
                (ax[i] - b[i]).abs() < 1e-8,
                "ax[{i}]={} b[{i}]={}",
                ax[i],
                b[i]
            );
        }
    }
    #[test]
    fn test_diagonal_preconditioner() {
        let a = CsrMatrix::from_triplets(3, 3, &[0, 1, 2], &[0, 1, 2], &[2.0, 4.0, 8.0]);
        let m = diagonal_preconditioner(&a);
        assert!((m[0] - 0.5).abs() < TOL);
        assert!((m[1] - 0.25).abs() < TOL);
        assert!((m[2] - 0.125).abs() < TOL);
    }
    #[test]
    fn test_diagonal_preconditioner_zero_diag() {
        let a = CsrMatrix::new(3, 3);
        let m = diagonal_preconditioner(&a);
        for v in m {
            assert!((v - 1.0).abs() < TOL);
        }
    }
    #[test]
    fn test_pcg_solve() {
        let lap = laplacian_1d(10, 1.0 / 11.0);
        let b = vec![1.0; 10];
        let precond = diagonal_preconditioner(&lap);
        let (x, _, res) = pcg_solve(&lap, &b, &precond, 1e-10, 500);
        let ax = lap.matvec(&x);
        let err: f64 = ax
            .iter()
            .zip(b.iter())
            .map(|(ai, bi)| (ai - bi).powi(2))
            .sum::<f64>()
            .sqrt();
        assert!(err < 1e-8, "residual={res}, err={err}");
    }
    #[test]
    fn test_nnz() {
        let m = laplacian_1d(5, 0.2);
        assert_eq!(m.nnz(), 5 + 2 * 4);
    }
    #[test]
    fn test_ilu0_solves_diagonal_system() {
        let a = CsrMatrix::from_triplets(3, 3, &[0, 1, 2], &[0, 1, 2], &[2.0, 4.0, 8.0]);
        let lu = SparseLu::ilu0(&a);
        let b = [4.0, 8.0, 16.0];
        let x = lu.solve(&b);
        let ax = a.matvec(&x);
        for i in 0..3 {
            assert!(
                (ax[i] - b[i]).abs() < 1e-8,
                "ILU0 solve: ax[{i}]={} b[{i}]={}",
                ax[i],
                b[i]
            );
        }
    }
    #[test]
    fn test_ilu0_spd_small() {
        let a = CsrMatrix::from_triplets(2, 2, &[0, 0, 1, 1], &[0, 1, 0, 1], &[4.0, 1.0, 1.0, 3.0]);
        let lu = SparseLu::ilu0(&a);
        let b = [1.0, 2.0];
        let x = lu.solve(&b);
        let ax = a.matvec(&x);
        for i in 0..2 {
            assert!(
                (ax[i] - b[i]).abs() < 1e-8,
                "ax[{i}]={} b[{i}]={}",
                ax[i],
                b[i]
            );
        }
    }
    #[test]
    fn test_ilu0_l_unit_diagonal() {
        let a = identity_csr(4);
        let lu = SparseLu::ilu0(&a);
        for i in 0..4 {
            assert!(
                (lu.l.get(i, i) - 1.0).abs() < 1e-12,
                "L diagonal should be 1"
            );
        }
    }
    #[test]
    fn test_ic0_spd_factor() {
        let a = CsrMatrix::from_triplets(
            3,
            3,
            &[0, 0, 1, 1, 1, 2, 2],
            &[0, 1, 0, 1, 2, 1, 2],
            &[4.0, -1.0, -1.0, 4.0, -1.0, -1.0, 4.0],
        );
        let ic = IncompleteCholesky::factor(&a).expect("IC should succeed for SPD matrix");
        for i in 0..3 {
            assert!(ic.l.get(i, i) > 0.0, "IC diagonal should be positive");
        }
    }
    #[test]
    fn test_ic0_apply_returns_vector() {
        let lap = laplacian_1d(5, 1.0 / 6.0);
        let ic = IncompleteCholesky::factor(&lap).expect("IC should succeed");
        let b = vec![1.0, 0.0, 0.0, 0.0, 0.0];
        let x = ic.apply(&b);
        assert_eq!(x.len(), 5);
        let nrm: f64 = x.iter().map(|v| v * v).sum::<f64>().sqrt();
        assert!(nrm > 0.0, "IC apply should return non-zero vector");
    }
    #[test]
    fn test_sparse_qr_least_squares_identity() {
        let id = identity_csr(3);
        let qr = SparseQr::factor(&id);
        let b = vec![1.0, 2.0, 3.0];
        let x = qr.solve_least_squares(&b);
        for i in 0..3 {
            assert!(
                (x[i] - b[i]).abs() < 1e-8,
                "QR solve: x[{i}]={} b[{i}]={}",
                x[i],
                b[i]
            );
        }
    }
    #[test]
    fn test_sparse_qr_overdetermined() {
        let a = CsrMatrix::from_triplets(3, 2, &[0, 1, 2, 2], &[0, 1, 0, 1], &[1.0, 1.0, 1.0, 1.0]);
        let qr = SparseQr::factor(&a);
        let b = vec![1.0, 1.0, 2.0];
        let x = qr.solve_least_squares(&b);
        assert_eq!(x.len(), 2);
        let ax = a.matvec(&x);
        let res: f64 = ax
            .iter()
            .zip(b.iter())
            .map(|(ai, bi)| (ai - bi).powi(2))
            .sum::<f64>();
        assert!(res < 1.0, "QR least-squares residual={res} should be small");
    }
    #[test]
    fn test_sparse_qr_q_orthogonal() {
        let a = laplacian_1d(3, 1.0);
        let qr = SparseQr::factor(&a);
        let m = qr.nrows;
        for i in 0..m {
            for j in 0..m {
                let dot: f64 = (0..m).map(|k| qr.q[i][k] * qr.q[j][k]).sum();
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (dot - expected).abs() < 1e-8,
                    "Q orthogonality failed at ({i},{j}): {dot}"
                );
            }
        }
    }
    #[test]
    fn test_ssor_apply_returns_correct_size() {
        let a = laplacian_1d(5, 1.0 / 6.0);
        let ssor = SsorPreconditioner::new(a, 1.0);
        let b = vec![1.0; 5];
        let x = ssor.apply(&b);
        assert_eq!(x.len(), 5, "SSOR apply should return vector of same size");
    }
    #[test]
    fn test_ssor_with_identity_omega1() {
        let a = identity_csr(4);
        let ssor = SsorPreconditioner::new(a, 1.0);
        let b = vec![2.0, 3.0, 4.0, 5.0];
        let x = ssor.apply(&b);
        assert_eq!(x.len(), 4);
        let nrm: f64 = x.iter().map(|v| v * v).sum::<f64>().sqrt();
        assert!(nrm > 0.0, "SSOR apply should return non-zero");
    }
    #[test]
    fn test_sparse_eig_power_identity() {
        let a = identity_csr(5);
        let result = sparse_eig_power(&a, 1e-10, 1000);
        let (lam, v) = result.expect("power iteration should converge");
        assert!(
            (lam - 1.0).abs() < 1e-6,
            "dominant eigenvalue of I should be 1, got {lam}"
        );
        let nrm: f64 = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        assert!(
            (nrm - 1.0).abs() < 1e-6,
            "eigenvector should be unit length"
        );
    }
    #[test]
    fn test_sparse_eig_power_diagonal() {
        let a = CsrMatrix::from_triplets(
            5,
            5,
            &[0, 1, 2, 3, 4],
            &[0, 1, 2, 3, 4],
            &[1.0, 2.0, 3.0, 4.0, 5.0],
        );
        let (lam, _) = sparse_eig_power(&a, 1e-10, 1000).expect("should converge");
        assert!(
            (lam - 5.0).abs() < 1e-6,
            "dominant eigenvalue should be 5, got {lam}"
        );
    }
    #[test]
    fn test_spectral_radius_estimate_identity() {
        let a = identity_csr(4);
        let sr = spectral_radius_estimate(&a, 1e-10, 1000);
        assert!(
            (sr - 1.0).abs() < 1e-4,
            "spectral radius of I should be 1, got {sr}"
        );
    }
    #[test]
    fn test_sparse_eig_smallest_laplacian() {
        let lap = laplacian_1d(5, 1.0 / 6.0);
        let result = sparse_eig_smallest(&lap, 1e-6, 500);
        if let Some((lam, _)) = result {
            assert!(
                lam > 0.0,
                "smallest eigenvalue of Laplacian should be positive, got {lam}"
            );
        }
    }
    #[test]
    fn test_csc_from_csr_identity() {
        let id = identity_csr(3);
        let csc = CscMatrix::from_csr(&id);
        assert_eq!(csc.nrows, 3);
        assert_eq!(csc.ncols, 3);
        for j in 0..3 {
            let col = csc.get_column(j);
            assert_eq!(col.len(), 1, "identity column {j} should have 1 entry");
            assert_eq!(col[0].0, j, "identity column {j} entry row should be {j}");
            assert!(
                (col[0].1 - 1.0).abs() < 1e-12,
                "identity column {j} value should be 1"
            );
        }
    }
    #[test]
    fn test_csc_column_extraction_laplacian() {
        let lap = laplacian_1d(5, 1.0);
        let csc = CscMatrix::from_csr(&lap);
        let col2 = csc.get_column(2);
        assert_eq!(
            col2.len(),
            3,
            "interior column of 1D Laplacian should have 3 entries"
        );
    }
    #[test]
    fn test_csc_column_extraction_values() {
        let a = CsrMatrix::from_triplets(3, 3, &[0, 0, 1, 2], &[0, 1, 1, 2], &[1.0, 2.0, 3.0, 4.0]);
        let csc = CscMatrix::from_csr(&a);
        let col1 = csc.get_column(1);
        assert_eq!(col1.len(), 2);
        let found_row0 = col1.iter().any(|&(r, v)| r == 0 && (v - 2.0).abs() < 1e-12);
        let found_row1 = col1.iter().any(|&(r, v)| r == 1 && (v - 3.0).abs() < 1e-12);
        assert!(
            found_row0,
            "column 1 should have row 0 entry with value 2.0"
        );
        assert!(
            found_row1,
            "column 1 should have row 1 entry with value 3.0"
        );
    }
    #[test]
    fn test_csc_matvec() {
        let id = identity_csr(4);
        let csc = CscMatrix::from_csr(&id);
        let x = vec![1.0, 2.0, 3.0, 4.0];
        let y = csc.matvec(&x);
        for i in 0..4 {
            assert!(
                (y[i] - x[i]).abs() < 1e-12,
                "CSC matvec identity: y[{i}]={} expected {}",
                y[i],
                x[i]
            );
        }
    }
    #[test]
    fn test_spgemm_identity_times_identity() {
        let id = identity_csr(3);
        let result = spgemm(&id, &id);
        let dense = result.to_dense();
        for (i, row) in dense.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (val - expected).abs() < 1e-12,
                    "I*I[{i}][{j}] = {} expected {}",
                    val,
                    expected
                );
            }
        }
    }
    #[test]
    fn test_spgemm_diagonal_times_diagonal() {
        let a = CsrMatrix::from_triplets(3, 3, &[0, 1, 2], &[0, 1, 2], &[2.0, 3.0, 4.0]);
        let b = CsrMatrix::from_triplets(3, 3, &[0, 1, 2], &[0, 1, 2], &[5.0, 6.0, 7.0]);
        let c = spgemm(&a, &b);
        let dense = c.to_dense();
        assert!(
            (dense[0][0] - 10.0).abs() < 1e-12,
            "C[0][0] = {} expected 10",
            dense[0][0]
        );
        assert!(
            (dense[1][1] - 18.0).abs() < 1e-12,
            "C[1][1] = {} expected 18",
            dense[1][1]
        );
        assert!(
            (dense[2][2] - 28.0).abs() < 1e-12,
            "C[2][2] = {} expected 28",
            dense[2][2]
        );
    }
    #[test]
    fn test_spgemm_rectangular() {
        let a = CsrMatrix::from_triplets(
            2,
            3,
            &[0, 0, 0, 1, 1, 1],
            &[0, 1, 2, 0, 1, 2],
            &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        );
        let b = CsrMatrix::from_triplets(
            3,
            2,
            &[0, 0, 1, 1, 2, 2],
            &[0, 1, 0, 1, 0, 1],
            &[7.0, 8.0, 9.0, 10.0, 11.0, 12.0],
        );
        let c = spgemm(&a, &b);
        assert_eq!(c.nrows, 2);
        assert_eq!(c.ncols, 2);
        let dense = c.to_dense();
        assert!(
            (dense[0][0] - 58.0).abs() < 1e-10,
            "C[0][0] = {} expected 58",
            dense[0][0]
        );
        assert!(
            (dense[0][1] - 64.0).abs() < 1e-10,
            "C[0][1] = {} expected 64",
            dense[0][1]
        );
        assert!(
            (dense[1][0] - 139.0).abs() < 1e-10,
            "C[1][0] = {} expected 139",
            dense[1][0]
        );
        assert!(
            (dense[1][1] - 154.0).abs() < 1e-10,
            "C[1][1] = {} expected 154",
            dense[1][1]
        );
    }
    #[test]
    fn test_spgemm_matches_matvec() {
        let lap = laplacian_1d(4, 1.0);
        let a2 = spgemm(&lap, &lap);
        let x = vec![1.0, 0.0, 0.0, 0.0];
        let direct = a2.matvec(&x);
        let two_step = lap.matvec(&lap.matvec(&x));
        for i in 0..4 {
            assert!(
                (direct[i] - two_step[i]).abs() < 1e-10,
                "SpGEMM result differs from sequential matvec at [{i}]"
            );
        }
    }
    #[test]
    fn test_coo_to_csr_basic() {
        let rows = vec![0usize, 1, 2];
        let cols = vec![0usize, 1, 2];
        let vals = vec![1.0, 2.0, 3.0];
        let csr = coo_to_csr(3, 3, &rows, &cols, &vals);
        assert_eq!(csr.nrows, 3);
        assert_eq!(csr.ncols, 3);
        assert!((csr.get(0, 0) - 1.0).abs() < 1e-12);
        assert!((csr.get(1, 1) - 2.0).abs() < 1e-12);
        assert!((csr.get(2, 2) - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_coo_to_csr_duplicate_sum() {
        let rows = vec![0usize, 0];
        let cols = vec![0usize, 0];
        let vals = vec![1.0, 2.0];
        let csr = coo_to_csr(2, 2, &rows, &cols, &vals);
        assert!(
            (csr.get(0, 0) - 3.0).abs() < 1e-12,
            "Duplicate entries should be summed"
        );
    }
    #[test]
    fn test_coo_to_csr_unsorted() {
        let rows = vec![2usize, 0, 1];
        let cols = vec![2usize, 0, 1];
        let vals = vec![3.0, 1.0, 2.0];
        let csr = coo_to_csr(3, 3, &rows, &cols, &vals);
        assert!(
            (csr.get(0, 0) - 1.0).abs() < 1e-12,
            "COO to CSR from unsorted: [0,0]"
        );
        assert!(
            (csr.get(1, 1) - 2.0).abs() < 1e-12,
            "COO to CSR from unsorted: [1,1]"
        );
        assert!(
            (csr.get(2, 2) - 3.0).abs() < 1e-12,
            "COO to CSR from unsorted: [2,2]"
        );
    }
    #[test]
    fn test_arnoldi_basis_orthonormal() {
        let a = laplacian_1d(6, 1.0);
        let b0 = vec![1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let (q, _h) = arnoldi(&a, &b0, 4);
        let k = q.len();
        for i in 0..k {
            for j in 0..k {
                let dot: f64 = q[i].iter().zip(q[j].iter()).map(|(a, b)| a * b).sum();
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (dot - expected).abs() < 1e-8,
                    "Arnoldi basis not orthonormal at ({i},{j}): dot = {dot}"
                );
            }
        }
    }
    #[test]
    fn test_arnoldi_krylov_relation() {
        let a = laplacian_1d(5, 1.0);
        let b0 = vec![1.0, 0.0, 0.0, 0.0, 0.0];
        let (q, h) = arnoldi(&a, &b0, 3);
        for j in 0..q.len().saturating_sub(1) {
            let aqj = a.matvec(&q[j]);
            let mut residual = aqj.clone();
            for i in 0..h.len() {
                if i < q.len() && j < h[i].len() {
                    for (r, qr) in residual.iter_mut().zip(q[i].iter()) {
                        *r -= h[i][j] * qr;
                    }
                }
            }
            let res_norm: f64 = residual.iter().map(|x| x * x).sum::<f64>().sqrt();
            assert!(
                res_norm < 1e-8,
                "Arnoldi relation violated at j={j}: residual norm = {res_norm}"
            );
        }
    }
    #[test]
    fn test_arnoldi_single_step_is_normalized_input() {
        let a = identity_csr(4);
        let b0 = vec![3.0, 4.0, 0.0, 0.0];
        let (q, _h) = arnoldi(&a, &b0, 1);
        assert_eq!(q.len(), 1);
        let expected = [0.6, 0.8, 0.0, 0.0];
        for i in 0..4 {
            assert!(
                (q[0][i] - expected[i]).abs() < 1e-10,
                "First Arnoldi vector: q[0][{i}] = {} expected {}",
                q[0][i],
                expected[i]
            );
        }
    }
    #[test]
    fn test_jacobi_scale_identity() {
        let id = identity_csr(4);
        let scaled = jacobi_scale(&id);
        for i in 0..4 {
            assert!(
                (scaled.get(i, i) - 1.0).abs() < 1e-12,
                "Jacobi-scaled identity diagonal[{i}] = {} expected 1",
                scaled.get(i, i)
            );
        }
    }
    #[test]
    fn test_jacobi_scale_diagonal_matrix() {
        let d = CsrMatrix::from_triplets(3, 3, &[0, 1, 2], &[0, 1, 2], &[2.0, 4.0, 8.0]);
        let scaled = jacobi_scale(&d);
        assert!((scaled.get(0, 0) - 1.0).abs() < 1e-12);
        assert!((scaled.get(1, 1) - 1.0).abs() < 1e-12);
        assert!((scaled.get(2, 2) - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_jacobi_scale_preserves_sparsity() {
        let lap = laplacian_1d(5, 1.0);
        let nnz_before = lap.nnz();
        let scaled = jacobi_scale(&lap);
        assert_eq!(
            scaled.nnz(),
            nnz_before,
            "Jacobi scaling should not change sparsity pattern"
        );
    }
    #[test]
    fn test_kronecker_identity_by_identity() {
        let i2 = identity_csr(2);
        let kron = i2.kronecker_product(&i2);
        assert_eq!(kron.nrows, 4);
        assert_eq!(kron.ncols, 4);
        for i in 0..4 {
            assert!(
                (kron.get(i, i) - 1.0).abs() < 1e-12,
                "diagonal mismatch at {}",
                i
            );
            for j in 0..4 {
                if i != j {
                    assert!(
                        kron.get(i, j).abs() < 1e-12,
                        "off-diagonal non-zero at ({},{})",
                        i,
                        j
                    );
                }
            }
        }
    }
    #[test]
    fn test_kronecker_product_dimensions() {
        let a = CsrMatrix::from_triplets(2, 3, &[0, 1], &[0, 2], &[1.0, 1.0]);
        let b = CsrMatrix::from_triplets(4, 5, &[0, 3], &[0, 4], &[1.0, 1.0]);
        let kron = a.kronecker_product(&b);
        assert_eq!(kron.nrows, 8);
        assert_eq!(kron.ncols, 15);
    }
    #[test]
    fn test_kronecker_scalar() {
        let a = CsrMatrix::from_triplets(1, 1, &[0], &[0], &[3.0]);
        let b = CsrMatrix::from_triplets(1, 1, &[0], &[0], &[4.0]);
        let kron = a.kronecker_product(&b);
        assert_eq!(kron.nrows, 1);
        assert_eq!(kron.ncols, 1);
        assert!((kron.get(0, 0) - 12.0).abs() < 1e-12);
    }
    #[test]
    fn test_ic0_identity_gives_identity() {
        let id = identity_csr(4);
        let l = id.incomplete_cholesky();
        for i in 0..4 {
            assert!(
                (l.get(i, i) - 1.0).abs() < 1e-10,
                "L[{},{}]={}",
                i,
                i,
                l.get(i, i)
            );
        }
    }
    #[test]
    fn test_ic0_diagonal_matrix() {
        let d = CsrMatrix::from_triplets(3, 3, &[0, 1, 2], &[0, 1, 2], &[4.0, 9.0, 16.0]);
        let l = d.incomplete_cholesky();
        assert!((l.get(0, 0) - 2.0).abs() < 1e-10);
        assert!((l.get(1, 1) - 3.0).abs() < 1e-10);
        assert!((l.get(2, 2) - 4.0).abs() < 1e-10);
    }
    #[test]
    fn test_ic0_lower_triangular() {
        let id = identity_csr(5);
        let l = id.incomplete_cholesky();
        for r in 0..5 {
            for c in (r + 1)..5 {
                assert!(
                    l.get(r, c).abs() < 1e-12,
                    "upper triangle non-zero L[{},{}]={}",
                    r,
                    c,
                    l.get(r, c)
                );
            }
        }
    }
    #[test]
    fn test_power_iteration_identity() {
        let id = identity_csr(4);
        let (lam, _v) = id.power_iteration(100, 1e-10);
        assert!((lam - 1.0).abs() < 1e-6, "lambda={}", lam);
    }
    #[test]
    fn test_power_iteration_diagonal() {
        let d =
            CsrMatrix::from_triplets(4, 4, &[0, 1, 2, 3], &[0, 1, 2, 3], &[1.0, 2.0, 3.0, 10.0]);
        let (lam, _v) = d.power_iteration(200, 1e-10);
        assert!((lam - 10.0).abs() < 1e-4, "expected lambda≈10, got {}", lam);
    }
    #[test]
    fn test_power_iteration_eigenvector_unit() {
        let id = identity_csr(3);
        let (_lam, v) = id.power_iteration(50, 1e-10);
        let norm: f64 = v.iter().map(|&x| x * x).sum::<f64>().sqrt();
        assert!(
            (norm - 1.0).abs() < 1e-9,
            "eigenvector not unit: norm={}",
            norm
        );
    }
}
