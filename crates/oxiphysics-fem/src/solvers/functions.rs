//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::sparse::CsrMatrix;

use super::types::{
    CholeskyDense, ConjugateGradient, ConvergenceMonitor, GaussSeidel, SolverStats,
};

/// Compute `y = A x` for a CSR sparse matrix.
///
/// `values`, `col_idx`, and `row_ptr` follow the standard CSR layout.
/// `x` has length `ncols`; the returned vector has length `nrows = row_ptr.len() - 1`.
pub fn sparse_matvec(values: &[f64], col_idx: &[usize], row_ptr: &[usize], x: &[f64]) -> Vec<f64> {
    let nrows = row_ptr.len().saturating_sub(1);
    let mut y = vec![0.0; nrows];
    for (i, y_i) in y.iter_mut().enumerate() {
        let mut sum = 0.0;
        for k in row_ptr[i]..row_ptr[i + 1] {
            sum += values[k] * x[col_idx[k]];
        }
        *y_i = sum;
    }
    y
}
/// Solve `A x = b` using Preconditioned Conjugate Gradient.
///
/// `precond[i] = 1 / A[i,i]` (Jacobi preconditioner).
///
/// Returns `(solution, stats)`.
pub fn pcg_solve(
    a: &CsrMatrix,
    b: &[f64],
    precond: &[f64],
    max_iter: usize,
    tol: f64,
) -> (Vec<f64>, SolverStats) {
    let n = b.len();
    let mut x = vec![0.0; n];
    let ax = a.mul_vec(&x);
    let mut r: Vec<f64> = b.iter().zip(ax.iter()).map(|(bi, ai)| bi - ai).collect();
    let mut z: Vec<f64> = r
        .iter()
        .zip(precond.iter())
        .map(|(ri, mi)| ri * mi)
        .collect();
    let mut p = z.clone();
    let mut rz_old = dot(&r, &z);
    let b_norm = dot(b, b).sqrt();
    let mut iters = 0;
    let mut converged = false;
    if b_norm < 1e-60 {
        let residual = dot(&r, &r).sqrt();
        return (
            x,
            SolverStats {
                iterations: 0,
                residual,
                converged: true,
            },
        );
    }
    for _ in 0..max_iter {
        iters += 1;
        let ap = a.mul_vec(&p);
        let p_ap = dot(&p, &ap);
        if p_ap.abs() < 1e-60 {
            break;
        }
        let alpha = rz_old / p_ap;
        for ((x_i, r_i), (&p_i, &ap_i)) in
            x.iter_mut().zip(r.iter_mut()).zip(p.iter().zip(ap.iter()))
        {
            *x_i += alpha * p_i;
            *r_i -= alpha * ap_i;
        }
        let r_norm = dot(&r, &r).sqrt();
        if r_norm / b_norm < tol {
            converged = true;
            break;
        }
        for ((z_i, &r_i), &pre_i) in z.iter_mut().zip(r.iter()).zip(precond.iter()) {
            *z_i = r_i * pre_i;
        }
        let rz_new = dot(&r, &z);
        let beta = rz_new / rz_old;
        for (p_i, &z_i) in p.iter_mut().zip(z.iter()) {
            *p_i = z_i + beta * *p_i;
        }
        rz_old = rz_new;
    }
    let residual = dot(&r, &r).sqrt();
    (
        x,
        SolverStats {
            iterations: iters,
            residual,
            converged,
        },
    )
}
/// Solve `A x = b` using the Conjugate Gradient method (sparse CSR).
pub fn conjugate_gradient(
    a: &CsrMatrix,
    b: &[f64],
    x0: &[f64],
    max_iter: usize,
    tol: f64,
) -> Vec<f64> {
    let n = b.len();
    assert_eq!(a.nrows, n);
    assert_eq!(a.ncols, n);
    assert_eq!(x0.len(), n);
    let mut x = x0.to_vec();
    let ax = a.mul_vec(&x);
    let mut r: Vec<f64> = b.iter().zip(ax.iter()).map(|(bi, ai)| bi - ai).collect();
    let mut p = r.clone();
    let mut rs_old = dot(&r, &r);
    if rs_old.sqrt() < tol {
        return x;
    }
    for _iter in 0..max_iter {
        let ap = a.mul_vec(&p);
        let p_ap = dot(&p, &ap);
        if p_ap.abs() < 1e-60 {
            break;
        }
        let alpha = rs_old / p_ap;
        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        let rs_new = dot(&r, &r);
        if rs_new.sqrt() < tol {
            break;
        }
        let beta = rs_new / rs_old;
        for i in 0..n {
            p[i] = r[i] + beta * p[i];
        }
        rs_old = rs_new;
    }
    x
}
/// Compute the Jacobi (diagonal) preconditioner for the given matrix.
pub fn jacobi_preconditioner(a: &CsrMatrix) -> Vec<f64> {
    let n = a.nrows;
    let mut precond = vec![1.0; n];
    for (i, precond_i) in precond.iter_mut().enumerate().take(n) {
        let diag = a.diagonal(i);
        if diag.abs() > 1e-30 {
            *precond_i = 1.0 / diag;
        }
    }
    precond
}
/// Solve `A x = b` using the Preconditioned Conjugate Gradient method (sparse CSR).
pub fn preconditioned_cg(
    a: &CsrMatrix,
    b: &[f64],
    precond: &[f64],
    max_iter: usize,
    tol: f64,
) -> Vec<f64> {
    let n = b.len();
    assert_eq!(a.nrows, n);
    assert_eq!(a.ncols, n);
    assert_eq!(precond.len(), n);
    let mut x = vec![0.0; n];
    let ax = a.mul_vec(&x);
    let mut r: Vec<f64> = b.iter().zip(ax.iter()).map(|(bi, ai)| bi - ai).collect();
    let mut z: Vec<f64> = r
        .iter()
        .zip(precond.iter())
        .map(|(ri, mi)| ri * mi)
        .collect();
    let mut p = z.clone();
    let mut rz_old = dot(&r, &z);
    let b_norm = dot(b, b).sqrt();
    if b_norm < 1e-60 {
        return x;
    }
    for _iter in 0..max_iter {
        let ap = a.mul_vec(&p);
        let p_ap = dot(&p, &ap);
        if p_ap.abs() < 1e-60 {
            break;
        }
        let alpha = rz_old / p_ap;
        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        let r_norm = dot(&r, &r).sqrt();
        if r_norm / b_norm < tol {
            break;
        }
        for i in 0..n {
            z[i] = r[i] * precond[i];
        }
        let rz_new = dot(&r, &z);
        let beta = rz_new / rz_old;
        for i in 0..n {
            p[i] = z[i] + beta * p[i];
        }
        rz_old = rz_new;
    }
    x
}
/// Dot product of two `CsrMatrix`-related slices.
pub(super) fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}
/// Dot product of plain slices.
pub(super) fn dot_slice(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}
/// Dense row-major matrix-vector product: y = A * x.
pub(super) fn dense_matvec(a: &[f64], x: &[f64], n: usize) -> Vec<f64> {
    let mut y = vec![0.0; n];
    for i in 0..n {
        for j in 0..n {
            y[i] += a[i * n + j] * x[j];
        }
    }
    y
}
/// Solve a dense square system `A x = b` (row-major, size `n`) using BiCGSTAB.
///
/// Returns `(solution, stats)`. Works for non-symmetric systems.
pub fn bicgstab_dense(
    a: &[f64],
    b: &[f64],
    n: usize,
    max_iter: usize,
    tol: f64,
) -> (Vec<f64>, SolverStats) {
    let mut x = vec![0.0_f64; n];
    let ax = dense_matvec(a, &x, n);
    let mut r: Vec<f64> = b.iter().zip(ax.iter()).map(|(bi, ai)| bi - ai).collect();
    let r_hat = r.clone();
    let b_norm = dot_slice(b, b).sqrt();
    let mut rho_old = 1.0_f64;
    let mut alpha = 1.0_f64;
    let mut omega = 1.0_f64;
    let mut p = vec![0.0_f64; n];
    let mut v = vec![0.0_f64; n];
    let mut iters = 0usize;
    let mut converged = false;
    if b_norm < 1e-60 {
        let residual = dot_slice(&r, &r).sqrt();
        return (
            x,
            SolverStats {
                iterations: 0,
                residual,
                converged: true,
            },
        );
    }
    for _ in 0..max_iter {
        iters += 1;
        let rho_new = dot_slice(&r_hat, &r);
        if rho_new.abs() < 1e-60 {
            break;
        }
        let beta = (rho_new / rho_old) * (alpha / omega);
        for i in 0..n {
            p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }
        v = dense_matvec(a, &p, n);
        let rv = dot_slice(&r_hat, &v);
        if rv.abs() < 1e-60 {
            break;
        }
        alpha = rho_new / rv;
        let mut s: Vec<f64> = r
            .iter()
            .zip(v.iter())
            .map(|(ri, vi)| ri - alpha * vi)
            .collect();
        let s_norm = dot_slice(&s, &s).sqrt();
        if s_norm / b_norm < tol {
            for i in 0..n {
                x[i] += alpha * p[i];
            }
            converged = true;
            r = s;
            break;
        }
        let t = dense_matvec(a, &s, n);
        let tt = dot_slice(&t, &t);
        omega = if tt.abs() > 1e-60 {
            dot_slice(&t, &s) / tt
        } else {
            0.0
        };
        for i in 0..n {
            x[i] += alpha * p[i] + omega * s[i];
            s[i] -= omega * t[i];
        }
        r = s;
        let r_norm = dot_slice(&r, &r).sqrt();
        if r_norm / b_norm < tol {
            converged = true;
            break;
        }
        rho_old = rho_new;
        if omega.abs() < 1e-60 {
            break;
        }
    }
    let residual = dot_slice(&r, &r).sqrt();
    (
        x,
        SolverStats {
            iterations: iters,
            residual,
            converged,
        },
    )
}
/// Solve a dense symmetric system `A x = b` using MINRES.
///
/// Suitable for indefinite symmetric systems (not just SPD).
/// Returns `(solution, stats)`.
pub fn minres_dense(
    a: &[f64],
    b: &[f64],
    n: usize,
    max_iter: usize,
    tol: f64,
) -> (Vec<f64>, SolverStats) {
    let mut x = vec![0.0_f64; n];
    let b_norm = dot_slice(b, b).sqrt();
    if b_norm < 1e-60 {
        return (
            x,
            SolverStats {
                iterations: 0,
                residual: 0.0,
                converged: true,
            },
        );
    }
    let mut v_prev = vec![0.0_f64; n];
    let mut v_curr: Vec<f64> = b.iter().map(|bi| bi / b_norm).collect();
    let mut beta_old = 0.0_f64;
    let mut beta_curr = b_norm;
    let mut c_old = 1.0_f64;
    let mut s_old = 0.0_f64;
    // c_cur and s_cur are always assigned before being read (first loop body),
    // so no meaningful initial value is needed.
    let mut c_cur;
    let mut s_cur;
    let mut eta = b_norm;
    let mut w = vec![0.0_f64; n];
    let mut w_old = vec![0.0_f64; n];
    let mut iters = 0usize;
    let mut converged = false;
    for _ in 0..max_iter {
        iters += 1;
        let av = dense_matvec(a, &v_curr, n);
        let alpha: f64 = dot_slice(&v_curr, &av);
        let mut v_next: Vec<f64> = av
            .iter()
            .zip(v_curr.iter())
            .zip(v_prev.iter())
            .map(|((avi, vci), vpi)| avi - alpha * vci - beta_curr * vpi)
            .collect();
        let beta_next = dot_slice(&v_next, &v_next).sqrt();
        if beta_next > 1e-60 {
            for vi in &mut v_next {
                *vi /= beta_next;
            }
        }
        let eps = s_old * beta_old;
        let delta = c_old * beta_old + s_old * alpha;
        let gamma_sq = (delta * delta + beta_next * beta_next).sqrt().max(1e-60);
        c_cur = delta / gamma_sq;
        s_cur = -beta_next / gamma_sq;
        for (i, w_i) in w.iter_mut().enumerate() {
            *w_i = (v_curr[i] - eps * w_old[i] - delta * *w_i) / gamma_sq;
        }
        for (x_i, &w_i) in x.iter_mut().zip(w.iter()) {
            *x_i += c_cur * eta * w_i;
        }
        eta *= -s_cur;
        let r_norm = eta.abs();
        if r_norm / b_norm < tol {
            converged = true;
            break;
        }
        v_prev = v_curr.clone();
        v_curr = v_next;
        beta_old = beta_curr;
        beta_curr = beta_next;
        c_old = c_cur;
        s_old = s_cur;
        w_old = w.clone();
        if beta_next < 1e-60 {
            converged = true;
            break;
        }
    }
    let ax = dense_matvec(a, &x, n);
    let r: Vec<f64> = b.iter().zip(ax.iter()).map(|(bi, ai)| bi - ai).collect();
    let residual = dot_slice(&r, &r).sqrt();
    if residual / b_norm < tol {
        converged = true;
    }
    (
        x,
        SolverStats {
            iterations: iters,
            residual,
            converged,
        },
    )
}
/// Solve a dense square system `A x = b` using GMRES(m) with restart.
///
/// `m` is the Krylov subspace dimension before restart.
/// Returns `(solution, stats)`.
pub fn gmres_dense(
    a: &[f64],
    b: &[f64],
    n: usize,
    m: usize,
    max_outer: usize,
    tol: f64,
) -> (Vec<f64>, SolverStats) {
    let mut x = vec![0.0_f64; n];
    let b_norm = dot_slice(b, b).sqrt();
    if b_norm < 1e-60 {
        return (
            x,
            SolverStats {
                iterations: 0,
                residual: 0.0,
                converged: true,
            },
        );
    }
    let m = m.max(1);
    let mut total_iters = 0usize;
    let mut converged = false;
    let mut final_residual = b_norm;
    for _ in 0..max_outer {
        let ax = dense_matvec(a, &x, n);
        let r0: Vec<f64> = b.iter().zip(ax.iter()).map(|(bi, ai)| bi - ai).collect();
        let r0_norm = dot_slice(&r0, &r0).sqrt();
        if r0_norm / b_norm < tol {
            converged = true;
            final_residual = r0_norm;
            break;
        }
        let mut v: Vec<Vec<f64>> = vec![vec![0.0; n]; m + 1];
        for i in 0..n {
            v[0][i] = r0[i] / r0_norm;
        }
        let mut h = vec![vec![0.0_f64; m]; m + 1];
        let mut g = vec![0.0_f64; m + 1];
        g[0] = r0_norm;
        let mut cs = vec![0.0_f64; m];
        let mut sn = vec![0.0_f64; m];
        let mut inner_iter = 0usize;
        for j in 0..m {
            inner_iter += 1;
            total_iters += 1;
            let mut w = dense_matvec(a, &v[j], n);
            for i in 0..=j {
                h[i][j] = dot_slice(&w, &v[i]);
                for k in 0..n {
                    w[k] -= h[i][j] * v[i][k];
                }
            }
            h[j + 1][j] = dot_slice(&w, &w).sqrt();
            if h[j + 1][j] > 1e-60 {
                for k in 0..n {
                    v[j + 1][k] = w[k] / h[j + 1][j];
                }
            }
            for i in 0..j {
                let tmp = cs[i] * h[i][j] - sn[i] * h[i + 1][j];
                h[i + 1][j] = sn[i] * h[i][j] + cs[i] * h[i + 1][j];
                h[i][j] = tmp;
            }
            let denom = (h[j][j] * h[j][j] + h[j + 1][j] * h[j + 1][j])
                .sqrt()
                .max(1e-60);
            cs[j] = h[j][j] / denom;
            sn[j] = -h[j + 1][j] / denom;
            h[j][j] = cs[j] * h[j][j] - sn[j] * h[j + 1][j];
            h[j + 1][j] = 0.0;
            g[j + 1] = sn[j] * g[j];
            g[j] *= cs[j];
            final_residual = g[j + 1].abs();
            if final_residual / b_norm < tol {
                converged = true;
                break;
            }
        }
        let sz = inner_iter;
        let mut y = vec![0.0_f64; sz];
        for i in (0..sz).rev() {
            y[i] = g[i];
            for k in (i + 1)..sz {
                y[i] -= h[i][k] * y[k];
            }
            if h[i][i].abs() > 1e-60 {
                y[i] /= h[i][i];
            }
        }
        for i in 0..sz {
            for k in 0..n {
                x[k] += y[i] * v[i][k];
            }
        }
        if converged {
            break;
        }
    }
    (
        x,
        SolverStats {
            iterations: total_iters,
            residual: final_residual,
            converged,
        },
    )
}
/// Compute incomplete LU factorization (ILU(0)) of a dense `n×n` matrix.
///
/// Returns `(L, U)` both as flat row-major vectors of size `n×n`.
/// L has unit diagonal; the sparsity pattern is preserved (ILU(0) = no fill-in).
pub fn ilu0_dense(a: &[f64], n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut l = vec![0.0_f64; n * n];
    let mut u = a.to_vec();
    for i in 0..n {
        l[i * n + i] = 1.0;
        for k in 0..i {
            if a[i * n + k].abs() < 1e-60 {
                continue;
            }
            if u[k * n + k].abs() < 1e-60 {
                continue;
            }
            let factor = u[i * n + k] / u[k * n + k];
            l[i * n + k] = factor;
            for j in k..n {
                u[i * n + j] -= factor * u[k * n + j];
            }
        }
    }
    (l, u)
}
/// Apply ILU(0) as a preconditioner: solve L y = b, then U x = y.
pub fn ilu0_solve(l: &[f64], u: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    let mut y = vec![0.0_f64; n];
    for i in 0..n {
        let mut s = b[i];
        for j in 0..i {
            s -= l[i * n + j] * y[j];
        }
        y[i] = s;
    }
    let mut x = vec![0.0_f64; n];
    for i in (0..n).rev() {
        let mut s = y[i];
        for j in (i + 1)..n {
            s -= u[i * n + j] * x[j];
        }
        x[i] = if u[i * n + i].abs() > 1e-60 {
            s / u[i * n + i]
        } else {
            0.0
        };
    }
    x
}
/// Compute the incomplete Cholesky factorization of a dense SPD `n×n` matrix.
///
/// Returns flat row-major lower-triangular `L` such that `L L^T ≈ A`.
/// Only entries where `A[i][j] != 0` are computed (ILU(0) pattern).
pub fn incomplete_cholesky_dense(a: &[f64], n: usize) -> Vec<f64> {
    let mut l = vec![0.0_f64; n * n];
    for i in 0..n {
        for j in 0..=i {
            if i != j && a[i * n + j].abs() < 1e-60 {
                continue;
            }
            let mut s = a[i * n + j];
            for k in 0..j {
                s -= l[i * n + k] * l[j * n + k];
            }
            if i == j {
                l[i * n + j] = if s > 0.0 { s.sqrt() } else { 1e-60_f64 };
            } else if l[j * n + j].abs() > 1e-60 {
                l[i * n + j] = s / l[j * n + j];
            }
        }
    }
    l
}
/// Compute the Schur complement `S = D - C A^{-1} B` for block systems.
///
/// `A`, `B`, `C`, `D` are all dense `n×n` row-major matrices.
/// Returns the Schur complement `S` as a flat `n×n` matrix.
pub fn schur_complement_dense(a: &[f64], b: &[f64], c: &[f64], d: &[f64], n: usize) -> Vec<f64> {
    let chol = CholeskyDense::new(n);
    let mut ainv_b = vec![0.0_f64; n * n];
    if let Ok(l) = chol.factorize(a) {
        for col in 0..n {
            let b_col: Vec<f64> = (0..n).map(|row| b[row * n + col]).collect();
            let x = chol.solve(&l, &b_col);
            for row in 0..n {
                ainv_b[row * n + col] = x[row];
            }
        }
    }
    let mut s = d.to_vec();
    for i in 0..n {
        for j in 0..n {
            let mut val = 0.0;
            for k in 0..n {
                val += c[i * n + k] * ainv_b[k * n + j];
            }
            s[i * n + j] -= val;
        }
    }
    s
}
/// Perform a two-level multigrid V-cycle for a dense system `A x = b`.
///
/// Uses Gauss-Seidel as smoother (`pre_smooth` and `post_smooth` iterations).
/// The coarse-grid operator is the lower-left quarter of `A`.
/// Returns `(solution, stats)`.
pub fn multigrid_vcycle_dense(
    a: &[f64],
    b: &[f64],
    _x0: &[f64],
    n: usize,
    pre_smooth: usize,
    post_smooth: usize,
) -> (Vec<f64>, SolverStats) {
    let gs = GaussSeidel::new(pre_smooth, 1e-60, 1.0);
    let mut x;
    x = gs.solve(a, b, n);
    let ax = dense_matvec(a, &x, n);
    let r: Vec<f64> = b.iter().zip(ax.iter()).map(|(bi, ai)| bi - ai).collect();
    let nc = (n / 2).max(1);
    let r_c: Vec<f64> = (0..nc).map(|i| r[i * 2]).collect();
    let mut a_c = vec![0.0_f64; nc * nc];
    for i in 0..nc {
        for j in 0..nc {
            a_c[i * nc + j] = a[(i * 2) * n + (j * 2)];
        }
    }
    let chol_c = CholeskyDense::new(nc);
    let e_c = if let Ok(lc) = chol_c.factorize(&a_c) {
        chol_c.solve(&lc, &r_c)
    } else {
        let cg = ConjugateGradient::new(nc * 10, 1e-10);
        cg.solve(&a_c, &r_c, nc)
    };
    for i in 0..nc {
        x[i * 2] += e_c[i];
        if i * 2 + 1 < n {
            x[i * 2 + 1] += 0.5 * (e_c[i] + if i + 1 < nc { e_c[i + 1] } else { 0.0 });
        }
    }
    let gs_post = GaussSeidel::new(post_smooth, 1e-60, 1.0);
    x = gs_post.solve(a, b, n);
    let ax_final = dense_matvec(a, &x, n);
    let r_final: Vec<f64> = b
        .iter()
        .zip(ax_final.iter())
        .map(|(bi, ai)| bi - ai)
        .collect();
    let residual = dot_slice(&r_final, &r_final).sqrt();
    let b_norm = dot_slice(b, b).sqrt();
    let converged = residual / b_norm.max(1e-60) < 1e-6;
    (
        x,
        SolverStats {
            iterations: pre_smooth + post_smooth,
            residual,
            converged,
        },
    )
}
/// Solve `A x = b` using Preconditioned CG with an ICC preconditioner.
///
/// The ICC preconditioner is computed from `a` automatically.  Returns
/// `(solution, stats)`.  Suitable for symmetric positive definite systems.
pub fn pcg_icc(a: &CsrMatrix, b: &[f64], max_iter: usize, tol: f64) -> (Vec<f64>, SolverStats) {
    use crate::sparse::IccPreconditioner;
    let icc = IccPreconditioner::new(a);
    let n = b.len();
    let mut x = vec![0.0f64; n];
    let mut r = b.to_vec();
    let mut z = icc.solve(&r);
    let mut p = z.clone();
    let mut rz_old = dot(&r, &z);
    let b_norm = dot(b, b).sqrt();
    let mut iters = 0usize;
    let mut converged = false;
    if b_norm < 1e-60 {
        let residual = dot(&r, &r).sqrt();
        return (
            x,
            SolverStats {
                iterations: 0,
                residual,
                converged: true,
            },
        );
    }
    for _ in 0..max_iter {
        iters += 1;
        let ap = a.mul_vec(&p);
        let p_ap = dot(&p, &ap);
        if p_ap.abs() < 1e-60 {
            break;
        }
        let alpha = rz_old / p_ap;
        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        let r_norm = dot(&r, &r).sqrt();
        if r_norm / b_norm < tol {
            converged = true;
            break;
        }
        z = icc.solve(&r);
        let rz_new = dot(&r, &z);
        let beta = rz_new / rz_old;
        for i in 0..n {
            p[i] = z[i] + beta * p[i];
        }
        rz_old = rz_new;
    }
    let residual = dot(&r, &r).sqrt();
    (
        x,
        SolverStats {
            iterations: iters,
            residual,
            converged,
        },
    )
}
/// Solve `A x = b` using PCG while recording the residual history in a
/// [`ConvergenceMonitor`].
///
/// Returns `(solution, monitor)`.
pub fn pcg_with_monitor(
    a: &CsrMatrix,
    b: &[f64],
    precond: &[f64],
    max_iter: usize,
    tol: f64,
) -> (Vec<f64>, ConvergenceMonitor) {
    let n = b.len();
    let mut monitor = ConvergenceMonitor::new(tol, max_iter);
    let mut x = vec![0.0f64; n];
    let mut r: Vec<f64> = b.to_vec();
    let mut z: Vec<f64> = r
        .iter()
        .zip(precond.iter())
        .map(|(ri, pi)| ri * pi)
        .collect();
    let mut p = z.clone();
    let mut rz_old = dot(&r, &z);
    let b_norm = dot(b, b).sqrt();
    if b_norm < 1e-60 {
        monitor.record(0.0);
        return (x, monitor);
    }
    for _ in 0..max_iter {
        let ap = a.mul_vec(&p);
        let p_ap = dot(&p, &ap);
        if p_ap.abs() < 1e-60 {
            break;
        }
        let alpha = rz_old / p_ap;
        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        let r_norm = dot(&r, &r).sqrt();
        monitor.record(r_norm / b_norm);
        if r_norm / b_norm < tol {
            break;
        }
        for i in 0..n {
            z[i] = r[i] * precond[i];
        }
        let rz_new = dot(&r, &z);
        let beta = rz_new / rz_old;
        for i in 0..n {
            p[i] = z[i] + beta * p[i];
        }
        rz_old = rz_new;
    }
    (x, monitor)
}
/// Solve a sparse system `A x = b` using BiCGSTAB (CSR format).
///
/// Suitable for non-symmetric or mildly indefinite systems.
/// Returns `(solution, stats)`.
pub fn bicgstab_sparse(
    a: &CsrMatrix,
    b: &[f64],
    max_iter: usize,
    tol: f64,
) -> (Vec<f64>, SolverStats) {
    let n = b.len();
    assert_eq!(a.nrows, n);
    assert_eq!(a.ncols, n);
    let mut x = vec![0.0f64; n];
    let ax = a.mul_vec(&x);
    let mut r: Vec<f64> = b.iter().zip(ax.iter()).map(|(bi, ai)| bi - ai).collect();
    let r_hat = r.clone();
    let b_norm = dot(b, b).sqrt();
    if b_norm < 1e-60 {
        let residual = dot(&r, &r).sqrt();
        return (
            x,
            SolverStats {
                iterations: 0,
                residual,
                converged: true,
            },
        );
    }
    let mut rho_old = 1.0f64;
    let mut alpha = 1.0f64;
    let mut omega = 1.0f64;
    let mut p = vec![0.0f64; n];
    let mut v = vec![0.0f64; n];
    let mut iters = 0usize;
    let mut converged = false;
    for _ in 0..max_iter {
        iters += 1;
        let rho_new = dot(&r_hat, &r);
        if rho_new.abs() < 1e-60 {
            break;
        }
        let beta = (rho_new / rho_old) * (alpha / omega);
        for i in 0..n {
            p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }
        v = a.mul_vec(&p);
        let rv = dot(&r_hat, &v);
        if rv.abs() < 1e-60 {
            break;
        }
        alpha = rho_new / rv;
        let mut s: Vec<f64> = r
            .iter()
            .zip(v.iter())
            .map(|(ri, vi)| ri - alpha * vi)
            .collect();
        let s_norm = dot(&s, &s).sqrt();
        if s_norm / b_norm < tol {
            for i in 0..n {
                x[i] += alpha * p[i];
            }
            converged = true;
            r = s;
            break;
        }
        let t = a.mul_vec(&s);
        let tt = dot(&t, &t);
        omega = if tt.abs() > 1e-60 {
            dot(&t, &s) / tt
        } else {
            0.0
        };
        for i in 0..n {
            x[i] += alpha * p[i] + omega * s[i];
            s[i] -= omega * t[i];
        }
        r = s;
        let r_norm = dot(&r, &r).sqrt();
        if r_norm / b_norm < tol {
            converged = true;
            break;
        }
        rho_old = rho_new;
        if omega.abs() < 1e-60 {
            break;
        }
    }
    let residual = dot(&r, &r).sqrt();
    (
        x,
        SolverStats {
            iterations: iters,
            residual,
            converged,
        },
    )
}
/// Solve a sparse symmetric system `A x = b` using MINRES.
///
/// Suitable for symmetric indefinite systems (SPD and SPSD are also handled).
/// Returns `(solution, stats)`.
pub fn minres_sparse(
    a: &CsrMatrix,
    b: &[f64],
    max_iter: usize,
    tol: f64,
) -> (Vec<f64>, SolverStats) {
    let n = b.len();
    assert_eq!(a.nrows, n);
    assert_eq!(a.ncols, n);
    let b_norm = dot(b, b).sqrt();
    let mut x = vec![0.0f64; n];
    if b_norm < 1e-60 {
        return (
            x,
            SolverStats {
                iterations: 0,
                residual: 0.0,
                converged: true,
            },
        );
    }
    let mut v_prev = vec![0.0f64; n];
    let mut v_curr: Vec<f64> = b.iter().map(|bi| bi / b_norm).collect();
    let mut beta_old = 0.0f64;
    let mut beta_curr = b_norm;
    let mut c_old = 1.0f64;
    let mut s_old = 0.0f64;
    // c_cur and s_cur are always assigned before being read (first loop body).
    let mut c_cur;
    let mut s_cur;
    let mut eta = b_norm;
    let mut w = vec![0.0f64; n];
    let mut w_old = vec![0.0f64; n];
    let mut iters = 0usize;
    let mut converged = false;
    for _ in 0..max_iter {
        iters += 1;
        let av = a.mul_vec(&v_curr);
        let alpha: f64 = dot(&v_curr, &av);
        let mut v_next: Vec<f64> = av
            .iter()
            .zip(v_curr.iter())
            .zip(v_prev.iter())
            .map(|((avi, vci), vpi)| avi - alpha * vci - beta_curr * vpi)
            .collect();
        let beta_next = dot(&v_next, &v_next).sqrt();
        if beta_next > 1e-60 {
            for vi in &mut v_next {
                *vi /= beta_next;
            }
        }
        let eps = s_old * beta_old;
        let delta = c_old * beta_old + s_old * alpha;
        let gamma_sq = (delta * delta + beta_next * beta_next).sqrt().max(1e-60);
        c_cur = delta / gamma_sq;
        s_cur = -beta_next / gamma_sq;
        for (i, w_i) in w.iter_mut().enumerate() {
            *w_i = (v_curr[i] - eps * w_old[i] - delta * *w_i) / gamma_sq;
        }
        for (x_i, &w_i) in x.iter_mut().zip(w.iter()) {
            *x_i += c_cur * eta * w_i;
        }
        eta *= -s_cur;
        let r_norm = eta.abs();
        if r_norm / b_norm < tol {
            converged = true;
            break;
        }
        v_prev = v_curr.clone();
        v_curr = v_next;
        beta_old = beta_curr;
        beta_curr = beta_next;
        c_old = c_cur;
        s_old = s_cur;
        w_old = w.clone();
        if beta_next < 1e-60 {
            converged = true;
            break;
        }
    }
    let ax_final = a.mul_vec(&x);
    let r: Vec<f64> = b
        .iter()
        .zip(ax_final.iter())
        .map(|(bi, ai)| bi - ai)
        .collect();
    let residual = dot(&r, &r).sqrt();
    if residual / b_norm < tol {
        converged = true;
    }
    (
        x,
        SolverStats {
            iterations: iters,
            residual,
            converged,
        },
    )
}
/// Solve a sparse system `A x = b` using GMRES(m) with restart.
///
/// `m` is the restart dimension (default recommendation: 20).
/// Returns `(solution, stats)`.
pub fn gmres_sparse(
    a: &CsrMatrix,
    b: &[f64],
    m: usize,
    max_outer: usize,
    tol: f64,
) -> (Vec<f64>, SolverStats) {
    let n = b.len();
    assert_eq!(a.nrows, n);
    assert_eq!(a.ncols, n);
    let b_norm = dot(b, b).sqrt();
    let mut x = vec![0.0f64; n];
    if b_norm < 1e-60 {
        return (
            x,
            SolverStats {
                iterations: 0,
                residual: 0.0,
                converged: true,
            },
        );
    }
    let m = m.max(1);
    let mut total_iters = 0usize;
    let mut converged = false;
    let mut final_residual = b_norm;
    for _ in 0..max_outer {
        let ax = a.mul_vec(&x);
        let r0: Vec<f64> = b.iter().zip(ax.iter()).map(|(bi, ai)| bi - ai).collect();
        let r0_norm = dot(&r0, &r0).sqrt();
        if r0_norm / b_norm < tol {
            converged = true;
            final_residual = r0_norm;
            break;
        }
        let mut v: Vec<Vec<f64>> = vec![vec![0.0; n]; m + 1];
        for i in 0..n {
            v[0][i] = r0[i] / r0_norm;
        }
        let mut h = vec![vec![0.0f64; m]; m + 1];
        let mut g = vec![0.0f64; m + 1];
        g[0] = r0_norm;
        let mut cs = vec![0.0f64; m];
        let mut sn = vec![0.0f64; m];
        let mut inner_iter = 0usize;
        for j in 0..m {
            inner_iter += 1;
            total_iters += 1;
            let mut w = a.mul_vec(&v[j]);
            for i in 0..=j {
                h[i][j] = dot(&w, &v[i]);
                let hi = h[i][j];
                for k in 0..n {
                    w[k] -= hi * v[i][k];
                }
            }
            h[j + 1][j] = dot(&w, &w).sqrt();
            if h[j + 1][j] > 1e-60 {
                let hj = h[j + 1][j];
                for k in 0..n {
                    v[j + 1][k] = w[k] / hj;
                }
            }
            for i in 0..j {
                let tmp = cs[i] * h[i][j] - sn[i] * h[i + 1][j];
                h[i + 1][j] = sn[i] * h[i][j] + cs[i] * h[i + 1][j];
                h[i][j] = tmp;
            }
            let denom = (h[j][j] * h[j][j] + h[j + 1][j] * h[j + 1][j])
                .sqrt()
                .max(1e-60);
            cs[j] = h[j][j] / denom;
            sn[j] = -h[j + 1][j] / denom;
            h[j][j] = cs[j] * h[j][j] - sn[j] * h[j + 1][j];
            h[j + 1][j] = 0.0;
            g[j + 1] = sn[j] * g[j];
            g[j] *= cs[j];
            final_residual = g[j + 1].abs();
            if final_residual / b_norm < tol {
                converged = true;
                break;
            }
        }
        let sz = inner_iter;
        let mut y = vec![0.0f64; sz];
        for i in (0..sz).rev() {
            y[i] = g[i];
            for k in (i + 1)..sz {
                y[i] -= h[i][k] * y[k];
            }
            if h[i][i].abs() > 1e-60 {
                y[i] /= h[i][i];
            }
        }
        for i in 0..sz {
            for k in 0..n {
                x[k] += y[i] * v[i][k];
            }
        }
        if converged {
            break;
        }
    }
    (
        x,
        SolverStats {
            iterations: total_iters,
            residual: final_residual,
            converged,
        },
    )
}
/// Solve `A x = b` using a simple two-level algebraic multigrid (AMG) V-cycle
/// as a preconditioner to PCG.
///
/// The coarse grid is built using the AMG coarsening/prolongation/Galerkin
/// operators from the `sparse` module.  This implements the outer PCG loop
/// with one AMG V-cycle per iteration.
///
/// Returns `(solution, stats)`.
pub fn amg_pcg(a: &CsrMatrix, b: &[f64], max_iter: usize, tol: f64) -> (Vec<f64>, SolverStats) {
    use crate::sparse::amg_v_cycle;
    let n = b.len();
    let mut x = vec![0.0f64; n];
    let ax = a.mul_vec(&x);
    let mut r: Vec<f64> = b.iter().zip(ax.iter()).map(|(bi, ai)| bi - ai).collect();
    let b_norm = dot(b, b).sqrt();
    if b_norm < 1e-60 {
        let residual = dot(&r, &r).sqrt();
        return (
            x,
            SolverStats {
                iterations: 0,
                residual,
                converged: true,
            },
        );
    }
    let mut z = amg_v_cycle(a, &r, &vec![0.0f64; n], 2, 2);
    let mut p = z.clone();
    let mut rz_old = dot(&r, &z);
    let mut iters = 0usize;
    let mut converged = false;
    for _ in 0..max_iter {
        iters += 1;
        let ap = a.mul_vec(&p);
        let p_ap = dot(&p, &ap);
        if p_ap.abs() < 1e-60 {
            break;
        }
        let alpha = rz_old / p_ap;
        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        let r_norm = dot(&r, &r).sqrt();
        if r_norm / b_norm < tol {
            converged = true;
            break;
        }
        z = amg_v_cycle(a, &r, &vec![0.0f64; n], 2, 2);
        let rz_new = dot(&r, &z);
        let beta = rz_new / rz_old;
        for i in 0..n {
            p[i] = z[i] + beta * p[i];
        }
        rz_old = rz_new;
    }
    let residual = dot(&r, &r).sqrt();
    (
        x,
        SolverStats {
            iterations: iters,
            residual,
            converged,
        },
    )
}
/// Nested iteration: solve `A x = b` starting from a coarse-grid approximation
/// then refine on the full system.
///
/// This performs `coarse_iter` iterations of PCG on a reduced system (first
/// `n_coarse` DOFs with a diagonal approximation), then uses the result as an
/// initial guess for `fine_iter` iterations of full PCG.
///
/// Returns `(solution, stats)`.
pub fn nested_iteration(
    a: &CsrMatrix,
    b: &[f64],
    n_coarse: usize,
    coarse_iter: usize,
    fine_iter: usize,
    tol: f64,
) -> (Vec<f64>, SolverStats) {
    let n = b.len();
    let n_c = n_coarse.min(n);
    let mut x = vec![0.0f64; n];
    if n_c > 0 {
        let mut coarse_triplets = Vec::new();
        for row in 0..n_c {
            let start = a.row_ptr[row];
            let end = a.row_ptr[row + 1];
            for idx in start..end {
                let col = a.col_indices[idx];
                if col < n_c {
                    coarse_triplets.push((row, col, a.values[idx]));
                }
            }
        }
        let a_c = CsrMatrix::from_triplets(n_c, n_c, &coarse_triplets);
        let b_c = b[..n_c].to_vec();
        let precond_c: Vec<f64> = (0..n_c)
            .map(|i| {
                let d = a_c.diagonal(i);
                if d.abs() > 1e-30 { 1.0 / d } else { 1.0 }
            })
            .collect();
        let (x_c, _) = pcg_solve(&a_c, &b_c, &precond_c, coarse_iter, tol * 10.0);
        x[..n_c].copy_from_slice(&x_c);
    }
    let precond: Vec<f64> = (0..n)
        .map(|i| {
            let d = a.diagonal(i);
            if d.abs() > 1e-30 { 1.0 / d } else { 1.0 }
        })
        .collect();
    let (x_fine, stats) = pcg_solve_from(a, b, &x, &precond, fine_iter, tol);
    (x_fine, stats)
}
/// PCG starting from a non-zero initial guess `x0`.
pub(super) fn pcg_solve_from(
    a: &CsrMatrix,
    b: &[f64],
    x0: &[f64],
    precond: &[f64],
    max_iter: usize,
    tol: f64,
) -> (Vec<f64>, SolverStats) {
    let n = b.len();
    let mut x = x0.to_vec();
    let ax = a.mul_vec(&x);
    let mut r: Vec<f64> = b.iter().zip(ax.iter()).map(|(bi, ai)| bi - ai).collect();
    let mut z: Vec<f64> = r
        .iter()
        .zip(precond.iter())
        .map(|(ri, pi)| ri * pi)
        .collect();
    let mut p = z.clone();
    let mut rz_old = dot(&r, &z);
    let b_norm = dot(b, b).sqrt();
    let mut iters = 0usize;
    let mut converged = false;
    if b_norm < 1e-60 {
        let residual = dot(&r, &r).sqrt();
        return (
            x,
            SolverStats {
                iterations: 0,
                residual,
                converged: true,
            },
        );
    }
    for _ in 0..max_iter {
        iters += 1;
        let ap = a.mul_vec(&p);
        let p_ap = dot(&p, &ap);
        if p_ap.abs() < 1e-60 {
            break;
        }
        let alpha = rz_old / p_ap;
        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        let r_norm = dot(&r, &r).sqrt();
        if r_norm / b_norm < tol {
            converged = true;
            break;
        }
        for i in 0..n {
            z[i] = r[i] * precond[i];
        }
        let rz_new = dot(&r, &z);
        let beta = rz_new / rz_old;
        for i in 0..n {
            p[i] = z[i] + beta * p[i];
        }
        rz_old = rz_new;
    }
    let residual = dot(&r, &r).sqrt();
    (
        x,
        SolverStats {
            iterations: iters,
            residual,
            converged,
        },
    )
}
/// Solve a dense square system `A x = b` via LU factorization with partial
/// pivoting.
///
/// Suitable for small systems (up to a few hundred DOFs).  For larger systems
/// prefer the iterative solvers.  Returns the solution vector, or a zero
/// vector if the system is singular.
pub fn lu_solve(a_in: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    assert_eq!(a_in.len(), n * n);
    assert_eq!(b.len(), n);
    let mut a = a_in.to_vec();
    let mut b = b.to_vec();
    let mut piv = vec![0usize; n];
    for col in 0..n {
        let mut max_val = a[col * n + col].abs();
        let mut max_row = col;
        for row in (col + 1)..n {
            let v = a[row * n + col].abs();
            if v > max_val {
                max_val = v;
                max_row = row;
            }
        }
        piv[col] = max_row;
        if max_row != col {
            for j in 0..n {
                a.swap(col * n + j, max_row * n + j);
            }
            b.swap(col, max_row);
        }
        let pivot = a[col * n + col];
        if pivot.abs() < 1e-60 {
            continue;
        }
        for row in (col + 1)..n {
            let factor = a[row * n + col] / pivot;
            a[row * n + col] = factor;
            for j in (col + 1)..n {
                let sub = factor * a[col * n + j];
                a[row * n + j] -= sub;
            }
            b[row] -= factor * b[col];
        }
    }
    let mut x = vec![0.0f64; n];
    for i in (0..n).rev() {
        let mut s = b[i];
        for j in (i + 1)..n {
            s -= a[i * n + j] * x[j];
        }
        x[i] = if a[i * n + i].abs() > 1e-60 {
            s / a[i * n + i]
        } else {
            0.0
        };
    }
    let _ = piv;
    x
}
