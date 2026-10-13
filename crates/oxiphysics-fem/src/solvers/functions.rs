//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::sparse::CsrMatrix;

use super::symmetric_amg::SymmetricAmg;
use super::types::{CholeskyDense, ConvergenceMonitor, GaussSeidel, SolverStats};

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
/// Suitable for indefinite symmetric systems (not just SPD).  The iteration is
/// the Paige–Saunders recurrence (Lanczos process plus a QR factorisation of
/// the tridiagonal matrix by Givens rotations), started from `x = 0`.
///
/// Returns `(solution, stats)`.  `stats.residual` is the true residual norm
/// `‖b − A x‖₂` of the returned vector and `stats.converged` is set only if
/// `‖b − A x‖₂ / ‖b‖₂ < tol`; it is never inferred from the recurrence.
pub fn minres_dense(
    a: &[f64],
    b: &[f64],
    n: usize,
    max_iter: usize,
    tol: f64,
) -> (Vec<f64>, SolverStats) {
    minres_iterate(b, max_iter, tol, |v| dense_matvec(a, v, n))
}
/// Shared MINRES iteration for a symmetric operator given as a mat-vec closure.
///
/// With Lanczos vectors `v_k`, coefficients `α_k`, `β_k` and the Givens
/// rotations `(c, s)` that triangularise the tridiagonal matrix, step `k`
/// computes
///
/// ```text
/// δ  = c₁ α − c₀ s₁ β_k           ρ₁ = √(δ² + β_{k+1}²)
/// ρ₂ = s₁ α + c₀ c₁ β_k           ρ₃ = s₀ β_k
/// w_k = (v_k − ρ₃ w_{k−2} − ρ₂ w_{k−1}) / ρ₁
/// x  += (δ/ρ₁) η w_k              η ← −(β_{k+1}/ρ₁) η
/// ```
///
/// where `(c₀, s₀)` and `(c₁, s₁)` are the rotations of steps `k−2` and `k−1`
/// and `|η|` is the norm of the current residual.  The three-term search
/// direction recurrence needs both `w_{k−1}` and `w_{k−2}`.
fn minres_iterate<F>(b: &[f64], max_iter: usize, tol: f64, matvec: F) -> (Vec<f64>, SolverStats)
where
    F: Fn(&[f64]) -> Vec<f64>,
{
    let n = b.len();
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
    let mut beta_curr = 0.0_f64;
    let (mut c_km2, mut s_km2) = (1.0_f64, 0.0_f64);
    let (mut c_km1, mut s_km1) = (1.0_f64, 0.0_f64);
    let mut eta = b_norm;
    let mut w_km2 = vec![0.0_f64; n];
    let mut w_km1 = vec![0.0_f64; n];
    let mut iters = 0usize;
    for _ in 0..max_iter {
        iters += 1;
        let av = matvec(&v_curr);
        let alpha = dot_slice(&v_curr, &av);
        let mut v_next: Vec<f64> = av
            .iter()
            .zip(v_curr.iter())
            .zip(v_prev.iter())
            .map(|((avi, vci), vpi)| avi - alpha * vci - beta_curr * vpi)
            .collect();
        let beta_next = dot_slice(&v_next, &v_next).sqrt();
        let delta = c_km1 * alpha - c_km2 * s_km1 * beta_curr;
        let rho1 = delta.hypot(beta_next);
        if rho1 <= 0.0 || rho1.is_nan() {
            break;
        }
        let rho2 = s_km1 * alpha + c_km2 * c_km1 * beta_curr;
        let rho3 = s_km2 * beta_curr;
        let c_new = delta / rho1;
        let s_new = beta_next / rho1;
        let w: Vec<f64> = v_curr
            .iter()
            .zip(w_km2.iter())
            .zip(w_km1.iter())
            .map(|((vi, w2), w1)| (vi - rho3 * w2 - rho2 * w1) / rho1)
            .collect();
        for (x_i, &w_i) in x.iter_mut().zip(w.iter()) {
            *x_i += c_new * eta * w_i;
        }
        eta *= -s_new;
        c_km2 = c_km1;
        s_km2 = s_km1;
        c_km1 = c_new;
        s_km1 = s_new;
        w_km2 = std::mem::replace(&mut w_km1, w);
        if eta.abs() / b_norm < tol {
            break;
        }
        // An invariant Krylov subspace was found: the iterate is already the
        // minimiser over the whole subspace and further steps add nothing.
        let av_norm = dot_slice(&av, &av).sqrt();
        if beta_next <= f64::EPSILON * av_norm {
            break;
        }
        for vi in &mut v_next {
            *vi /= beta_next;
        }
        v_prev = std::mem::replace(&mut v_curr, v_next);
        beta_curr = beta_next;
    }
    let ax = matvec(&x);
    let residual = b
        .iter()
        .zip(ax.iter())
        .map(|(bi, ai)| (bi - ai) * (bi - ai))
        .sum::<f64>()
        .sqrt();
    let converged = residual / b_norm < tol;
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
/// Linear-interpolation prolongation for a line of `n` unknowns with
/// full coarsening: coarse node `i` coincides with fine node `2i + 1`.
///
/// Returns, for every fine node, its `(coarse index, weight)` entries: a coarse
/// node is injected with weight 1, and a fine node between two coarse nodes
/// takes the weight ½ from each neighbour (a missing neighbour beyond a
/// homogeneous Dirichlet end contributes nothing).
fn linear_interpolation_rows(n: usize, n_coarse: usize) -> Vec<Vec<(usize, f64)>> {
    (0..n)
        .map(|j| {
            if j % 2 == 1 {
                vec![((j - 1) / 2, 1.0)]
            } else {
                let mut entries = Vec::with_capacity(2);
                if j >= 2 {
                    entries.push((j / 2 - 1, 0.5));
                }
                if j / 2 < n_coarse {
                    entries.push((j / 2, 0.5));
                }
                entries
            }
        })
        .collect()
}
/// Perform one two-grid V-cycle for a dense system `A x = b`.
///
/// The cycle starts from the initial guess `x0` (the zero vector if
/// `x0.len() != n`) and consists of
///
/// 1. `pre_smooth` Gauss-Seidel sweeps applied **to the current iterate**,
/// 2. a coarse-grid correction with the Galerkin operator `A_c = Pᵀ A P`,
///    where `P` is linear interpolation from every second unknown (the
///    unknowns are assumed to be ordered along a line, as for a 1-D mesh),
///    solved exactly by Cholesky (pivoted LU if `A_c` is not SPD), and
/// 3. `post_smooth` Gauss-Seidel sweeps applied to the corrected iterate.
///
/// Returns `(solution, stats)`.  `stats.converged` is `true` when the final
/// relative residual `‖b − A x‖₂ / ‖b‖₂` is below `1e-6`.
pub fn multigrid_vcycle_dense(
    a: &[f64],
    b: &[f64],
    x0: &[f64],
    n: usize,
    pre_smooth: usize,
    post_smooth: usize,
) -> (Vec<f64>, SolverStats) {
    let mut x = if x0.len() == n {
        x0.to_vec()
    } else {
        vec![0.0_f64; n]
    };
    GaussSeidel::new(pre_smooth, 1e-60, 1.0).smooth_in_place(a, b, n, &mut x);
    let nc = n / 2;
    if nc > 0 {
        let ax = dense_matvec(a, &x, n);
        let r: Vec<f64> = b.iter().zip(ax.iter()).map(|(bi, ai)| bi - ai).collect();
        let p = linear_interpolation_rows(n, nc);
        let mut r_c = vec![0.0_f64; nc];
        for (r_j, row) in r.iter().zip(p.iter()) {
            for &(ci, w) in row {
                r_c[ci] += w * r_j;
            }
        }
        let mut a_c = vec![0.0_f64; nc * nc];
        for (j, row_j) in p.iter().enumerate() {
            for (k, row_k) in p.iter().enumerate() {
                let a_jk = a[j * n + k];
                if a_jk == 0.0 {
                    continue;
                }
                for &(ci, wi) in row_j {
                    for &(cj, wj) in row_k {
                        a_c[ci * nc + cj] += wi * a_jk * wj;
                    }
                }
            }
        }
        let chol_c = CholeskyDense::new(nc);
        let e_c = match chol_c.factorize(&a_c) {
            Ok(lc) => chol_c.solve(&lc, &r_c),
            Err(_) => lu_solve(&a_c, &r_c, nc),
        };
        for (x_j, row) in x.iter_mut().zip(p.iter()) {
            for &(ci, w) in row {
                *x_j += w * e_c[ci];
            }
        }
    }
    GaussSeidel::new(post_smooth, 1e-60, 1.0).smooth_in_place(a, b, n, &mut x);
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
/// See [`minres_dense`] for the recurrence and the meaning of the returned
/// statistics: `stats.residual` is the true residual and `stats.converged`
/// reflects only `‖b − A x‖₂ / ‖b‖₂ < tol`.
///
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
    minres_iterate(b, max_iter, tol, |v| a.mul_vec(v))
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
/// Solve `A x = b` with PCG preconditioned by an algebraic multigrid V-cycle.
///
/// The grid hierarchy (coarsening, interpolation and Galerkin coarse
/// operators) is built **once** by [`SymmetricAmg`] before the iteration
/// starts, and every PCG iteration applies one symmetric positive definite
/// V-cycle (identical damped-Jacobi pre- and post-smoothing, exact coarsest
/// solve), as the conjugate-gradient method requires of its preconditioner.
///
/// If the hierarchy cannot be built (for example because `A` is only
/// positive semi-definite, so the coarsest Cholesky factorisation breaks down)
/// the solver falls back to Jacobi-preconditioned CG.
///
/// Returns `(solution, stats)`.
pub fn amg_pcg(a: &CsrMatrix, b: &[f64], max_iter: usize, tol: f64) -> (Vec<f64>, SolverStats) {
    let n = b.len();
    assert_eq!(a.nrows, n);
    assert_eq!(a.ncols, n);
    let mut x = vec![0.0f64; n];
    let mut r: Vec<f64> = b.to_vec();
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
    let amg = match SymmetricAmg::new(a) {
        Ok(amg) => amg,
        Err(_) => return pcg_solve(a, b, &jacobi_preconditioner(a), max_iter, tol),
    };
    let mut z = vec![0.0f64; n];
    amg.apply(&r, &mut z);
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
        amg.apply(&r, &mut z);
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
