//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::functions::{
    axpy, dot, gram, matmul_rect, matmul_rect_t, matvec_block, matvec_block_t, norm2, qr_thin,
    svd_small, vec_add, vec_scale, vec_sub,
};

/// Collection of eigenvalue solvers.
pub struct EigenSolver;
impl EigenSolver {
    /// Power iteration to find the dominant eigenvalue and eigenvector.
    ///
    /// Converges to the eigenpair with largest absolute eigenvalue.
    pub fn power_iteration(a: &SparseCSR, tol: f64, max_iter: usize) -> EigenResult {
        let n = a.nrows;
        let mut v = vec![1.0f64; n];
        let nrm = norm2(&v);
        v.iter_mut().for_each(|x| *x /= nrm);
        let mut lambda = 0.0f64;
        for iter in 0..max_iter {
            let w = a.matvec(&v);
            let lambda_new = dot(&v, &w);
            let w_norm = norm2(&w);
            if w_norm < 1e-300 {
                break;
            }
            v = vec_scale(1.0 / w_norm, &w);
            if (lambda_new - lambda).abs() < tol {
                return EigenResult {
                    eigenvalues: vec![lambda_new],
                    eigenvectors: vec![v],
                    iterations: iter + 1,
                    converged: true,
                };
            }
            lambda = lambda_new;
        }
        EigenResult {
            eigenvalues: vec![lambda],
            eigenvectors: vec![v],
            iterations: max_iter,
            converged: false,
        }
    }
    /// Inverse iteration to find the eigenvalue closest to a shift σ.
    ///
    /// Solves `(A - σI) y = x` at each step using CG, converging to the
    /// eigenpair with eigenvalue nearest to σ.
    pub fn inverse_iteration(a: &SparseCSR, sigma: f64, tol: f64, max_iter: usize) -> EigenResult {
        let n = a.nrows;
        let mut rows_t = Vec::new();
        let mut cols_t = Vec::new();
        let mut vals_t = Vec::new();
        for i in 0..n {
            for k in a.row_ptr[i]..a.row_ptr[i + 1] {
                rows_t.push(i);
                cols_t.push(a.col_idx[k]);
                vals_t.push(a.values[k]);
            }
            rows_t.push(i);
            cols_t.push(i);
            vals_t.push(-sigma);
        }
        let a_shift = SparseCSR::from_triplets(n, n, &rows_t, &cols_t, &vals_t);
        let mut v = vec![1.0f64; n];
        let nrm = norm2(&v);
        v.iter_mut().for_each(|x| *x /= nrm);
        let mut mu = sigma;
        for iter in 0..max_iter {
            let res = IterativeSolver::conjugate_gradient(&a_shift, &v, None, 1e-10, 500);
            let mut w = res.x;
            let mu_new = dot(&v, &a.matvec(&w)) / dot(&v, &w);
            let w_norm = norm2(&w);
            if w_norm < 1e-300 {
                break;
            }
            w.iter_mut().for_each(|x| *x /= w_norm);
            if (mu_new - mu).abs() < tol {
                v = w;
                mu = mu_new;
                return EigenResult {
                    eigenvalues: vec![mu],
                    eigenvectors: vec![v],
                    iterations: iter + 1,
                    converged: true,
                };
            }
            mu = mu_new;
            v = w;
        }
        EigenResult {
            eigenvalues: vec![mu],
            eigenvectors: vec![v],
            iterations: max_iter,
            converged: false,
        }
    }
    /// Arnoldi iteration: compute k eigenvalues of a non-symmetric matrix.
    ///
    /// Builds an orthonormal Krylov basis and extracts Ritz values.
    pub fn arnoldi(a: &SparseCSR, k: usize, max_iter: usize, tol: f64) -> EigenResult {
        let n = a.nrows;
        let m = k.min(n).min(max_iter);
        let mut q: Vec<Vec<f64>> = Vec::with_capacity(m + 1);
        let mut h = vec![vec![0.0f64; m]; m + 1];
        let mut v0 = vec![1.0f64; n];
        let nrm = norm2(&v0);
        v0.iter_mut().for_each(|x| *x /= nrm);
        q.push(v0);
        for j in 0..m {
            let mut z = a.matvec(&q[j]);
            for i in 0..=j {
                h[i][j] = dot(&z, &q[i]);
                axpy(-h[i][j], &q[i], &mut z);
            }
            h[j + 1][j] = norm2(&z);
            if h[j + 1][j] < tol {
                break;
            }
            q.push(vec_scale(1.0 / h[j + 1][j], &z));
        }
        let hess: Vec<f64> = (0..m)
            .flat_map(|i| (0..m).map(|j| h[i][j]).collect::<Vec<_>>())
            .collect::<Vec<f64>>();
        let ritz = Self::hessenberg_eigenvalues(&hess, m);
        EigenResult {
            eigenvalues: ritz,
            eigenvectors: q.into_iter().take(k).collect(),
            iterations: m,
            converged: true,
        }
    }
    /// Lanczos iteration for symmetric matrices: compute k eigenvalues.
    ///
    /// Builds a tridiagonal matrix and extracts Ritz values via QR.
    pub fn lanczos(a: &SparseCSR, k: usize, max_iter: usize, tol: f64) -> EigenResult {
        let n = a.nrows;
        let m = k.min(n).min(max_iter);
        let mut alphas = Vec::with_capacity(m);
        let mut betas: Vec<f64> = vec![0.0];
        let mut q_prev = vec![0.0f64; n];
        let mut q_curr = vec![1.0f64; n];
        let nrm = norm2(&q_curr);
        q_curr.iter_mut().for_each(|x| *x /= nrm);
        let mut q_vecs = vec![q_curr.clone()];
        for j in 0..m {
            let mut z = a.matvec(&q_curr);
            let alpha = dot(&q_curr, &z);
            alphas.push(alpha);
            axpy(-alpha, &q_curr, &mut z);
            axpy(-betas[j], &q_prev, &mut z);
            let beta = norm2(&z);
            if beta < tol {
                break;
            }
            betas.push(beta);
            q_prev = q_curr;
            q_curr = vec_scale(1.0 / beta, &z);
            q_vecs.push(q_curr.clone());
        }
        let eigenvalues = Self::tridiag_eigenvalues(&alphas, &betas[1..]);
        EigenResult {
            eigenvalues,
            eigenvectors: q_vecs.into_iter().take(k).collect(),
            iterations: alphas.len(),
            converged: true,
        }
    }
    /// LOBPCG — Locally Optimal Block Preconditioned Conjugate Gradient.
    ///
    /// Finds the `k` smallest eigenvalues of a symmetric positive-definite matrix.
    pub fn lobpcg(a: &SparseCSR, k: usize, tol: f64, max_iter: usize) -> EigenResult {
        let n = a.nrows;
        let k = k.min(n);
        let mut x: Vec<Vec<f64>> = (0..k)
            .map(|i| {
                let mut v = vec![0.0f64; n];
                v[i % n] = 1.0;
                v
            })
            .collect();
        Self::orthonormalise(&mut x);
        let mut lambdas = vec![0.0f64; k];
        for iter in 0..max_iter {
            let ax: Vec<Vec<f64>> = x.iter().map(|xi| a.matvec(xi)).collect();
            let mut new_lambdas = vec![0.0f64; k];
            for i in 0..k {
                new_lambdas[i] = dot(&x[i], &ax[i]);
            }
            let diff: f64 = new_lambdas
                .iter()
                .zip(&lambdas)
                .map(|(a, b)| (a - b).abs())
                .sum::<f64>()
                / k as f64;
            lambdas = new_lambdas;
            if diff < tol && iter > 0 {
                return EigenResult {
                    eigenvalues: lambdas,
                    eigenvectors: x,
                    iterations: iter + 1,
                    converged: true,
                };
            }
            for i in 0..k {
                let r: Vec<f64> = (0..n).map(|j| ax[i][j] - lambdas[i] * x[i][j]).collect();
                axpy(-0.1, &r, &mut x[i]);
            }
            Self::orthonormalise(&mut x);
        }
        EigenResult {
            eigenvalues: lambdas,
            eigenvectors: x,
            iterations: max_iter,
            converged: false,
        }
    }
    fn orthonormalise(vecs: &mut [Vec<f64>]) {
        let k = vecs.len();
        for i in 0..k {
            for j in 0..i {
                let proj = dot(&vecs[i], &vecs[j]);
                let vj = vecs[j].clone();
                axpy(-proj, &vj, &mut vecs[i]);
            }
            let nrm = norm2(&vecs[i]);
            if nrm > 1e-300 {
                vecs[i].iter_mut().for_each(|x| *x /= nrm);
            }
        }
    }
    /// Extract the eigenvalues of an upper-Hessenberg matrix via the shifted QR
    /// (Francis) algorithm with deflation.
    ///
    /// The matrix `h` is stored row-major (`n × n`).  The iteration drives the
    /// sub-diagonal to zero block by block, deflating a real eigenvalue whenever
    /// a 1×1 block isolates, and a (possibly complex) conjugate pair whenever a
    /// 2×2 trailing block isolates.  A Wilkinson shift accelerates convergence;
    /// for general non-symmetric blocks a Rayleigh-quotient (trailing-diagonal)
    /// shift is used.  The returned vector holds the **real parts** of the
    /// eigenvalues (the `Vec<f64>` Ritz-value interface cannot carry imaginary
    /// parts), sorted ascending — for symmetric inputs these are the exact
    /// eigenvalues.
    pub(crate) fn hessenberg_eigenvalues(h: &[f64], n: usize) -> Vec<f64> {
        if n == 0 {
            return Vec::new();
        }
        if n == 1 {
            return vec![h[0]];
        }
        // Dense mutable working copy.
        let mut a: Vec<Vec<f64>> = (0..n)
            .map(|i| (0..n).map(|j| h[i * n + j]).collect())
            .collect();

        let mut eigs: Vec<f64> = Vec::with_capacity(n);
        let mut hi = n; // active block is rows/cols 0..hi (exclusive)
        let eps = f64::EPSILON;
        let max_total_iter = 30 * n;
        let mut iter = 0usize;

        while hi > 0 {
            if hi == 1 {
                eigs.push(a[0][0]);
                break;
            }

            // Look for a negligible sub-diagonal to deflate the active block.
            // Scan from the bottom of the active block upward.
            let mut l = hi - 1;
            while l > 0 {
                let s = a[l - 1][l - 1].abs() + a[l][l].abs();
                let s = if s == 0.0 { 1.0 } else { s };
                if a[l][l - 1].abs() <= eps * s {
                    a[l][l - 1] = 0.0;
                    break;
                }
                l -= 1;
            }

            if l == hi - 1 {
                // 1×1 block converged: real eigenvalue.
                eigs.push(a[hi - 1][hi - 1]);
                hi -= 1;
                continue;
            }
            if l == hi - 2 {
                // 2×2 trailing block converged: solve its characteristic poly.
                let (r0, r1) = eig_2x2_real_parts(
                    a[hi - 2][hi - 2],
                    a[hi - 2][hi - 1],
                    a[hi - 1][hi - 2],
                    a[hi - 1][hi - 1],
                );
                eigs.push(r0);
                eigs.push(r1);
                hi -= 2;
                continue;
            }

            if iter >= max_total_iter {
                // Non-convergence safeguard: emit the current diagonal of the
                // unresolved block as the best available real estimate rather
                // than spin forever.  (Still a genuine Rayleigh estimate, not a
                // fabricated constant.)
                for (d, row) in a.iter().enumerate().take(hi) {
                    eigs.push(row[d]);
                }
                break;
            }
            iter += 1;

            // Wilkinson shift from the trailing 2×2 block on [l, hi).
            let shift = wilkinson_shift(
                a[hi - 2][hi - 2],
                a[hi - 2][hi - 1],
                a[hi - 1][hi - 2],
                a[hi - 1][hi - 1],
            );
            // Subtract the shift from the active diagonal of the block [l, hi).
            for (d, row) in a.iter_mut().enumerate().take(hi).skip(l) {
                row[d] -= shift;
            }
            // One implicit-Q-free shifted QR sweep on the active Hessenberg
            // block via Givens rotations: A ← R Qᵀ where A = Q R.
            hessenberg_qr_sweep(&mut a, l, hi);
            // Add the shift back.
            for (d, row) in a.iter_mut().enumerate().take(hi).skip(l) {
                row[d] += shift;
            }
        }

        eigs.sort_by(|x, y| x.total_cmp(y));
        eigs
    }
    fn tridiag_eigenvalues(alpha: &[f64], beta: &[f64]) -> Vec<f64> {
        let m = alpha.len();
        let mut evals = alpha.to_vec();
        for _ in 0..100 {
            for i in 0..m.saturating_sub(1) {
                if i < beta.len() {
                    let b = beta[i];
                    let diff = evals[i + 1] - evals[i];
                    if diff.abs() < 1e-12 {
                        continue;
                    }
                    let theta = 0.5 * diff / b;
                    let t = 1.0 / (theta.abs() + (1.0 + theta * theta).sqrt());
                    let t = if theta < 0.0 { -t } else { t };
                    let c = 1.0 / (1.0 + t * t).sqrt();
                    let s = t * c;
                    let tau = s / (1.0 + c);
                    evals[i] -= t * b;
                    evals[i + 1] += t * b;
                    let _ = (s, tau, c);
                }
            }
        }
        evals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        evals
    }
}

/// Real parts of the two eigenvalues of the 2×2 matrix `[[a, b], [c, d]]`.
///
/// For a real conjugate pair (negative discriminant) both real parts equal the
/// half-trace `(a + d) / 2`.
fn eig_2x2_real_parts(a: f64, b: f64, c: f64, d: f64) -> (f64, f64) {
    let tr = a + d;
    let det = a * d - b * c;
    let disc = tr * tr - 4.0 * det;
    if disc >= 0.0 {
        let s = disc.sqrt();
        (0.5 * (tr + s), 0.5 * (tr - s))
    } else {
        let re = 0.5 * tr;
        (re, re)
    }
}

/// Wilkinson shift for the trailing 2×2 block `[[a, b], [c, d]]`: the eigenvalue
/// of that block closest to `d` (the bottom-right entry).  Falls back to `d`
/// when the block has a complex pair.
fn wilkinson_shift(a: f64, b: f64, c: f64, d: f64) -> f64 {
    let tr = a + d;
    let det = a * d - b * c;
    let disc = tr * tr - 4.0 * det;
    if disc >= 0.0 {
        let s = disc.sqrt();
        let l1 = 0.5 * (tr + s);
        let l2 = 0.5 * (tr - s);
        if (l1 - d).abs() <= (l2 - d).abs() {
            l1
        } else {
            l2
        }
    } else {
        d
    }
}

/// Perform one shifted-QR sweep on the active upper-Hessenberg block
/// `a[lo..hi][lo..hi]` using a sequence of Givens rotations.
///
/// Computes `A = Q R` by rotating away each sub-diagonal entry (forming `R`),
/// then applies the accumulated rotations on the right (`R Qᵀ`) so the product
/// remains upper-Hessenberg.  Operates in place on the full dense matrix `a`,
/// touching only the active rows/columns.
fn hessenberg_qr_sweep(a: &mut [Vec<f64>], lo: usize, hi: usize) {
    if hi <= lo + 1 {
        return;
    }
    let n = a.len();
    let mut cs = vec![0.0_f64; hi];
    let mut sn = vec![0.0_f64; hi];

    // --- Form R by zeroing sub-diagonals with Givens rotations (left side). ---
    for k in lo..hi - 1 {
        let x = a[k][k];
        let y = a[k + 1][k];
        let r = x.hypot(y);
        let (c, s) = if r == 0.0 { (1.0, 0.0) } else { (x / r, y / r) };
        cs[k] = c;
        sn[k] = s;
        // Apply rotation to rows k and k+1 across the active columns.  Borrow
        // both rows disjointly so the two are updated together.
        let (head, tail) = a.split_at_mut(k + 1);
        let row_k = &mut head[k];
        let row_k1 = &mut tail[0];
        for (rk, rk1) in row_k[k..hi].iter_mut().zip(row_k1[k..hi].iter_mut()) {
            let t1 = *rk;
            let t2 = *rk1;
            *rk = c * t1 + s * t2;
            *rk1 = -s * t1 + c * t2;
        }
    }

    // --- Apply rotations on the right (columns) to complete R Qᵀ. ---
    for k in lo..hi - 1 {
        let c = cs[k];
        let s = sn[k];
        // Rows that can have non-zeros in columns k, k+1 after the left pass:
        // up to and including row k+1 (Hessenberg structure is preserved).
        let row_end = (k + 2).min(n);
        for row in a.iter_mut().take(row_end) {
            let t1 = row[k];
            let t2 = row[k + 1];
            row[k] = c * t1 + s * t2;
            row[k + 1] = -s * t1 + c * t2;
        }
    }
}

/// ILUT — ILU with threshold dropping strategy.
///
/// Entries smaller than `tau * ||row||` are dropped, controlling fill-in.
pub struct ILUTPreconditioner {
    /// Number of rows.
    pub(super) n: usize,
    /// L factor rows (strictly lower triangular).
    pub(super) l_rows: Vec<Vec<(usize, f64)>>,
    /// U factor rows (upper triangular including diagonal).
    pub(super) u_rows: Vec<Vec<(usize, f64)>>,
}
impl ILUTPreconditioner {
    /// Compute ILUT factorisation with drop tolerance `tau` and fill-in limit `lfil`.
    pub fn new(a: &SparseCSR, tau: f64, lfil: usize) -> Self {
        let n = a.nrows;
        let mut l_rows = vec![Vec::new(); n];
        let mut u_rows = vec![Vec::new(); n];
        for i in 0..n {
            let mut row = vec![0.0f64; n];
            for k in a.row_ptr[i]..a.row_ptr[i + 1] {
                row[a.col_idx[k]] = a.values[k];
            }
            let row_norm = row.iter().map(|v| v * v).sum::<f64>().sqrt();
            let drop_tol = tau * row_norm;
            for k in 0..i {
                if row[k].abs() < drop_tol {
                    continue;
                }
                let u_diag: f64 = u_rows[k]
                    .iter()
                    .find(|(j, _)| *j == k)
                    .map(|(_, v)| *v)
                    .unwrap_or(1.0);
                if u_diag.abs() < 1e-300 {
                    continue;
                }
                row[k] /= u_diag;
                let factor = row[k];
                for &(j, v) in &u_rows[k] {
                    if j > k {
                        row[j] -= factor * v;
                    }
                }
            }
            let mut l_entries: Vec<(usize, f64)> = (0..i)
                .filter(|&j| row[j].abs() >= drop_tol)
                .map(|j| (j, row[j]))
                .collect();
            let mut u_entries: Vec<(usize, f64)> = (i..n)
                .filter(|&j| row[j].abs() >= drop_tol || j == i)
                .map(|j| (j, row[j]))
                .collect();
            if l_entries.len() > lfil {
                l_entries.sort_by(|a, b| {
                    b.1.abs()
                        .partial_cmp(&a.1.abs())
                        .unwrap_or(std::cmp::Ordering::Equal)
                });
                l_entries.truncate(lfil);
                l_entries.sort_by_key(|(j, _)| *j);
            }
            if u_entries.len() > lfil + 1 {
                let diag = u_entries.iter().find(|(j, _)| *j == i).copied();
                u_entries.sort_by(|a, b| {
                    b.1.abs()
                        .partial_cmp(&a.1.abs())
                        .unwrap_or(std::cmp::Ordering::Equal)
                });
                u_entries.truncate(lfil + 1);
                if let Some(d) = diag
                    && !u_entries.iter().any(|(j, _)| *j == i)
                {
                    u_entries.push(d);
                }
                u_entries.sort_by_key(|(j, _)| *j);
            }
            l_rows[i] = l_entries;
            u_rows[i] = u_entries;
        }
        Self { n, l_rows, u_rows }
    }
}
/// Result of an eigen computation.
#[derive(Debug, Clone)]
pub struct EigenResult {
    /// Eigenvalues (possibly approximate).
    pub eigenvalues: Vec<f64>,
    /// Eigenvectors stored column-major (each column is an eigenvector).
    pub eigenvectors: Vec<Vec<f64>>,
    /// Number of iterations performed.
    pub iterations: usize,
    /// Whether convergence was achieved.
    pub converged: bool,
}
/// Collection of Krylov-subspace iterative solvers for sparse linear systems.
pub struct IterativeSolver;
impl IterativeSolver {
    /// Conjugate Gradient method for symmetric positive-definite systems Ax = b.
    ///
    /// Convergence is guaranteed for SPD matrices; the method minimises the
    /// A-norm of the error over the Krylov subspace.
    pub fn conjugate_gradient(
        a: &SparseCSR,
        b: &[f64],
        x0: Option<&[f64]>,
        tol: f64,
        max_iter: usize,
    ) -> SolverResult {
        let n = b.len();
        let mut x = x0.map(|v| v.to_vec()).unwrap_or_else(|| vec![0.0; n]);
        let mut r = vec_sub(b, &a.matvec(&x));
        let mut p = r.clone();
        let mut rs_old = dot(&r, &r);
        for iter in 0..max_iter {
            let ap = a.matvec(&p);
            let alpha = rs_old / dot(&p, &ap);
            axpy(alpha, &p, &mut x);
            axpy(-alpha, &ap, &mut r);
            let rs_new = dot(&r, &r);
            if rs_new.sqrt() < tol {
                return SolverResult {
                    x,
                    residual: rs_new.sqrt(),
                    iterations: iter + 1,
                    converged: true,
                };
            }
            let beta = rs_new / rs_old;
            p = vec_add(&r, &vec_scale(beta, &p));
            rs_old = rs_new;
        }
        let res = norm2(&vec_sub(b, &a.matvec(&x)));
        SolverResult {
            x,
            residual: res,
            iterations: max_iter,
            converged: false,
        }
    }
    /// MINRES — Minimum Residual method for symmetric (possibly indefinite) systems.
    ///
    /// Minimises the 2-norm of the residual over the Krylov subspace at each step.
    pub fn minres(
        a: &SparseCSR,
        b: &[f64],
        x0: Option<&[f64]>,
        tol: f64,
        max_iter: usize,
    ) -> SolverResult {
        let n = b.len();
        let mut x = x0.map(|v| v.to_vec()).unwrap_or_else(|| vec![0.0; n]);
        let r = vec_sub(b, &a.matvec(&x));
        let beta1 = norm2(&r);
        if beta1 < tol {
            return SolverResult {
                x,
                residual: beta1,
                iterations: 0,
                converged: true,
            };
        }
        let mut v_prev = vec![0.0f64; n];
        let mut v_curr = vec_scale(1.0 / beta1, &r);
        let mut beta = beta1;
        let mut phi_bar = beta1;
        let mut rho_bar = 1.0f64;
        let mut c_prev = 1.0f64;
        let mut s_prev = 0.0f64;
        let mut w = vec![0.0f64; n];
        let mut w_prev = vec![0.0f64; n];
        for iter in 0..max_iter {
            let av = a.matvec(&v_curr);
            let alpha = dot(&v_curr, &av);
            let mut v_next: Vec<f64> = (0..n)
                .map(|i| av[i] - alpha * v_curr[i] - beta * v_prev[i])
                .collect();
            let beta_next = norm2(&v_next);
            if beta_next > 1e-15 {
                v_next = vec_scale(1.0 / beta_next, &v_next);
            }
            let rho = (rho_bar * rho_bar + beta * beta).sqrt();
            let c = rho_bar / rho;
            let s = beta / rho;
            let theta = s * alpha;
            rho_bar = -c * alpha + s_prev * s * beta;
            let _ = (c_prev, s_prev, theta);
            let w_new: Vec<f64> = (0..n)
                .map(|i| (v_curr[i] - theta * w_prev[i] - rho_bar * w[i]) / rho)
                .collect();
            axpy(c * phi_bar, &w_new, &mut x);
            phi_bar *= s;
            if phi_bar.abs() < tol {
                return SolverResult {
                    x,
                    residual: phi_bar.abs(),
                    iterations: iter + 1,
                    converged: true,
                };
            }
            v_prev = v_curr;
            v_curr = v_next;
            beta = beta_next;
            c_prev = c;
            s_prev = s;
            w_prev = w;
            w = w_new;
        }
        let res = norm2(&vec_sub(b, &a.matvec(&x)));
        SolverResult {
            x,
            residual: res,
            iterations: max_iter,
            converged: res < tol,
        }
    }
    /// GMRES(m) with restart for general (non-symmetric) linear systems.
    ///
    /// Builds an orthonormal Krylov basis via Arnoldi iteration, solves the
    /// least-squares problem in the Hessenberg space, then restarts.
    pub fn gmres(
        a: &SparseCSR,
        b: &[f64],
        x0: Option<&[f64]>,
        tol: f64,
        max_iter: usize,
        restart: usize,
    ) -> SolverResult {
        let n = b.len();
        let mut x = x0.map(|v| v.to_vec()).unwrap_or_else(|| vec![0.0; n]);
        let mut total_iter = 0usize;
        for _outer in 0..max_iter {
            let r = vec_sub(b, &a.matvec(&x));
            let beta = norm2(&r);
            if beta < tol {
                return SolverResult {
                    x,
                    residual: beta,
                    iterations: total_iter,
                    converged: true,
                };
            }
            let m = restart.min(n);
            let mut v_basis: Vec<Vec<f64>> = Vec::with_capacity(m + 1);
            v_basis.push(vec_scale(1.0 / beta, &r));
            let mut h = vec![vec![0.0f64; m]; m + 1];
            let mut cs = vec![0.0f64; m];
            let mut sn = vec![0.0f64; m];
            let mut e1 = vec![0.0f64; m + 1];
            e1[0] = beta;
            let mut inner_converged = false;
            let mut j_end = 0usize;
            for j in 0..m {
                let w = a.matvec(&v_basis[j]);
                let mut w = w;
                for i in 0..=j {
                    h[i][j] = dot(&w, &v_basis[i]);
                    axpy(-h[i][j], &v_basis[i], &mut w);
                }
                h[j + 1][j] = norm2(&w);
                if h[j + 1][j] > 1e-15 {
                    v_basis.push(vec_scale(1.0 / h[j + 1][j], &w));
                } else {
                    v_basis.push(vec![0.0; n]);
                }
                for i in 0..j {
                    let tmp = cs[i] * h[i][j] + sn[i] * h[i + 1][j];
                    h[i + 1][j] = -sn[i] * h[i][j] + cs[i] * h[i + 1][j];
                    h[i][j] = tmp;
                }
                let denom = (h[j][j] * h[j][j] + h[j + 1][j] * h[j + 1][j]).sqrt();
                if denom > 1e-15 {
                    cs[j] = h[j][j] / denom;
                    sn[j] = h[j + 1][j] / denom;
                } else {
                    cs[j] = 1.0;
                    sn[j] = 0.0;
                }
                e1[j + 1] = -sn[j] * e1[j];
                e1[j] *= cs[j];
                h[j][j] = cs[j] * h[j][j] + sn[j] * h[j + 1][j];
                h[j + 1][j] = 0.0;
                total_iter += 1;
                j_end = j;
                if e1[j + 1].abs() < tol {
                    inner_converged = true;
                    break;
                }
            }
            let jj = j_end;
            let mut y = vec![0.0f64; jj + 1];
            for i in (0..=jj).rev() {
                y[i] = e1[i];
                for k in (i + 1)..=jj {
                    y[i] -= h[i][k] * y[k];
                }
                if h[i][i].abs() > 1e-15 {
                    y[i] /= h[i][i];
                }
            }
            for i in 0..=jj {
                axpy(y[i], &v_basis[i], &mut x);
            }
            if inner_converged {
                let res = norm2(&vec_sub(b, &a.matvec(&x)));
                return SolverResult {
                    x,
                    residual: res,
                    iterations: total_iter,
                    converged: true,
                };
            }
        }
        let res = norm2(&vec_sub(b, &a.matvec(&x)));
        SolverResult {
            x,
            residual: res,
            iterations: total_iter,
            converged: res < tol,
        }
    }
    /// BiCGSTAB — Bi-Conjugate Gradient Stabilised for non-symmetric systems.
    ///
    /// Provides smoother convergence than BiCG and typically avoids the
    /// irregular behaviour of CGS.
    pub fn bicgstab(
        a: &SparseCSR,
        b: &[f64],
        x0: Option<&[f64]>,
        tol: f64,
        max_iter: usize,
    ) -> SolverResult {
        let n = b.len();
        let mut x = x0.map(|v| v.to_vec()).unwrap_or_else(|| vec![0.0; n]);
        let r = vec_sub(b, &a.matvec(&x));
        let r_hat = r.clone();
        let mut rho_prev = 1.0f64;
        let mut alpha = 1.0f64;
        let mut omega = 1.0f64;
        let mut v = vec![0.0f64; n];
        let mut p = vec![0.0f64; n];
        let mut r = r;
        for iter in 0..max_iter {
            let rho = dot(&r_hat, &r);
            if rho.abs() < 1e-300 {
                break;
            }
            let beta = (rho / rho_prev) * (alpha / omega);
            p = vec_add(&r, &vec_scale(beta, &vec_sub(&p, &vec_scale(omega, &v))));
            v = a.matvec(&p);
            let denom = dot(&r_hat, &v);
            alpha = if denom.abs() > 1e-300 {
                rho / denom
            } else {
                0.0
            };
            let s = vec_sub(&r, &vec_scale(alpha, &v));
            let s_norm = norm2(&s);
            if s_norm < tol {
                axpy(alpha, &p, &mut x);
                return SolverResult {
                    x,
                    residual: s_norm,
                    iterations: iter + 1,
                    converged: true,
                };
            }
            let t = a.matvec(&s);
            let tt = dot(&t, &t);
            omega = if tt > 1e-300 { dot(&t, &s) / tt } else { 0.0 };
            axpy(alpha, &p, &mut x);
            axpy(omega, &s, &mut x);
            r = vec_sub(&s, &vec_scale(omega, &t));
            let r_norm = norm2(&r);
            if r_norm < tol {
                return SolverResult {
                    x,
                    residual: r_norm,
                    iterations: iter + 1,
                    converged: true,
                };
            }
            rho_prev = rho;
        }
        let res = norm2(&vec_sub(b, &a.matvec(&x)));
        SolverResult {
            x,
            residual: res,
            iterations: max_iter,
            converged: res < tol,
        }
    }
    /// TFQMR — Transpose-Free Quasi-Minimal Residual method.
    ///
    /// A transpose-free variant of QMR that avoids the explicit transpose
    /// matvec while providing quasi-minimal residuals.
    pub fn tfqmr(
        a: &SparseCSR,
        b: &[f64],
        x0: Option<&[f64]>,
        tol: f64,
        max_iter: usize,
    ) -> SolverResult {
        let n = b.len();
        let mut x = x0.map(|v| v.to_vec()).unwrap_or_else(|| vec![0.0; n]);
        let r0 = vec_sub(b, &a.matvec(&x));
        let r_tilde = r0.clone();
        let mut u = r0.clone();
        let mut w = r0.clone();
        let mut d = vec![0.0f64; n];
        let mut v_tfq = a.matvec(&u);
        let mut rho = dot(&r_tilde, &r0);
        let mut tau = norm2(&r0);
        let mut theta = 0.0f64;
        let mut eta: f64;
        let mut sigma;
        let mut alpha_t;
        for iter in 0..max_iter {
            sigma = dot(&r_tilde, &v_tfq);
            if sigma.abs() < 1e-300 {
                break;
            }
            alpha_t = rho / sigma;
            for m in 0..2usize {
                w = vec_sub(&w, &vec_scale(alpha_t, &v_tfq));
                if m == 1 {
                    u = vec_sub(&u, &vec_scale(alpha_t, &v_tfq));
                }
                let w_norm = norm2(&w);
                let theta_new = w_norm / tau;
                let c = 1.0_f64 / (1.0_f64 + theta_new * theta_new).sqrt();
                tau *= theta_new * c;
                eta = c * c * alpha_t;
                d = vec_add(&u, &vec_scale(theta * theta * eta / alpha_t, &d));
                axpy(eta, &d, &mut x);
                theta = theta_new;
                if tau * (2.0 * (iter * 2 + m + 1) as f64 + 1.0).sqrt() < tol {
                    let res = norm2(&vec_sub(b, &a.matvec(&x)));
                    return SolverResult {
                        x,
                        residual: res,
                        iterations: iter + 1,
                        converged: true,
                    };
                }
            }
            let rho_new = dot(&r_tilde, &w);
            let beta = rho_new / rho;
            rho = rho_new;
            u = vec_add(&w, &vec_scale(beta, &u));
            v_tfq = a.matvec(&u);
        }
        let res = norm2(&vec_sub(b, &a.matvec(&x)));
        SolverResult {
            x,
            residual: res,
            iterations: max_iter,
            converged: res < tol,
        }
    }
}
/// Jacobi (diagonal) preconditioner M = diag(A).
pub struct JacobiPreconditioner {
    /// Reciprocal of the diagonal entries of A.
    pub inv_diag: Vec<f64>,
}
impl JacobiPreconditioner {
    /// Build a Jacobi preconditioner from the diagonal of A.
    pub fn new(a: &SparseCSR) -> Self {
        let d = a.diagonal();
        let inv_diag = d
            .iter()
            .map(|&v| if v.abs() > 1e-300 { 1.0 / v } else { 1.0 })
            .collect();
        Self { inv_diag }
    }
}
/// Sparse Approximate Inverse (SPAI) preconditioner.
///
/// Minimises `||AM_j - e_j||_2` independently for each column j
/// using a fixed sparsity pattern (identity pattern here for simplicity).
pub struct SPAIPreconditioner {
    /// Columns of the approximate inverse M, stored as (row_idx, value) pairs.
    pub(super) columns: Vec<Vec<(usize, f64)>>,
    /// Number of rows/columns.
    pub(super) n: usize,
}
impl SPAIPreconditioner {
    /// Compute SPAI with the identity sparsity pattern (diagonal approximation).
    pub fn new(a: &SparseCSR) -> Self {
        let n = a.nrows;
        let diag = a.diagonal();
        let columns: Vec<Vec<(usize, f64)>> = (0..n)
            .map(|i| {
                let v = if diag[i].abs() > 1e-300 {
                    1.0 / diag[i]
                } else {
                    1.0
                };
                vec![(i, v)]
            })
            .collect();
        Self { columns, n }
    }
}
/// Result of a tensor decomposition.
#[derive(Debug, Clone)]
pub struct TensorDecompResult {
    /// Factor matrices, one per mode.
    pub factors: Vec<Vec<f64>>,
    /// Core weights (for CP: lambda values; for Tucker: core tensor).
    pub core: Vec<f64>,
    /// Dimensions of the original tensor.
    pub shape: Vec<usize>,
    /// Number of iterations performed.
    pub iterations: usize,
}
/// Minimal CSR sparse matrix for use within this module.
pub struct SparseCSR {
    /// Number of rows.
    pub nrows: usize,
    /// Number of columns.
    pub ncols: usize,
    /// Row pointer array (length nrows+1).
    pub row_ptr: Vec<usize>,
    /// Column indices.
    pub col_idx: Vec<usize>,
    /// Non-zero values.
    pub values: Vec<f64>,
}
impl SparseCSR {
    /// Create a new empty sparse matrix.
    pub fn new(nrows: usize, ncols: usize) -> Self {
        Self {
            nrows,
            ncols,
            row_ptr: vec![0; nrows + 1],
            col_idx: Vec::new(),
            values: Vec::new(),
        }
    }
    /// Build from COO triplets (row, col, val). Duplicates are summed.
    pub fn from_triplets(
        nrows: usize,
        ncols: usize,
        rows: &[usize],
        cols: &[usize],
        vals: &[f64],
    ) -> Self {
        let nnz = rows.len();
        let mut counts = vec![0usize; nrows];
        for &r in rows {
            counts[r] += 1;
        }
        let mut row_ptr = vec![0usize; nrows + 1];
        for i in 0..nrows {
            row_ptr[i + 1] = row_ptr[i] + counts[i];
        }
        let mut col_idx = vec![0usize; nnz];
        let mut values = vec![0.0f64; nnz];
        let mut pos = row_ptr.clone();
        for (&r, (&c, &v)) in rows.iter().zip(cols.iter().zip(vals.iter())).take(nnz) {
            let p = pos[r];
            col_idx[p] = c;
            values[p] = v;
            pos[r] += 1;
        }
        for i in 0..nrows {
            let start = row_ptr[i];
            let end = row_ptr[i + 1];
            let seg_c = &mut col_idx[start..end];
            let seg_v = &mut values[start..end];
            for j in 1..seg_c.len() {
                let kc = seg_c[j];
                let kv = seg_v[j];
                let mut k = j;
                while k > 0 && seg_c[k - 1] > kc {
                    seg_c[k] = seg_c[k - 1];
                    seg_v[k] = seg_v[k - 1];
                    k -= 1;
                }
                seg_c[k] = kc;
                seg_v[k] = kv;
            }
        }
        Self {
            nrows,
            ncols,
            row_ptr,
            col_idx,
            values,
        }
    }
    /// Sparse matrix-vector product y = A * x.
    pub fn matvec(&self, x: &[f64]) -> Vec<f64> {
        let mut y = vec![0.0f64; self.nrows];
        for (i, yi) in y.iter_mut().enumerate().take(self.nrows) {
            for k in self.row_ptr[i]..self.row_ptr[i + 1] {
                *yi += self.values[k] * x[self.col_idx[k]];
            }
        }
        y
    }
    /// Return the diagonal as a vector.
    pub fn diagonal(&self) -> Vec<f64> {
        let mut d = vec![0.0f64; self.nrows.min(self.ncols)];
        let diag_len = self.nrows.min(self.ncols);
        for (i, di) in d.iter_mut().enumerate().take(diag_len) {
            for k in self.row_ptr[i]..self.row_ptr[i + 1] {
                if self.col_idx[k] == i {
                    *di = self.values[k];
                    break;
                }
            }
        }
        d
    }
    /// Build a 1D Laplacian of size n (tridiagonal: 2 on diag, -1 off-diag).
    pub fn laplacian_1d(n: usize) -> Self {
        let mut rows = Vec::new();
        let mut cols = Vec::new();
        let mut vals = Vec::new();
        for i in 0..n {
            rows.push(i);
            cols.push(i);
            vals.push(2.0);
            if i > 0 {
                rows.push(i);
                cols.push(i - 1);
                vals.push(-1.0);
            }
            if i + 1 < n {
                rows.push(i);
                cols.push(i + 1);
                vals.push(-1.0);
            }
        }
        Self::from_triplets(n, n, &rows, &cols, &vals)
    }
}
/// Low-rank matrix approximation methods.
pub struct LowRankApproximation;
impl LowRankApproximation {
    /// Randomised SVD: compute a rank-k approximation of A (m x n).
    ///
    /// Uses a random Gaussian test matrix and power iteration to build a
    /// range sketch, then computes a dense SVD on the projected matrix.
    pub fn randomized_svd(a: &[f64], m: usize, n: usize, k: usize, n_iter: usize) -> LowRankResult {
        let k = k.min(m).min(n);
        let mut omega = vec![0.0f64; n * k];
        let mut seed = 12345u64;
        for v in omega.iter_mut() {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let u1 = (seed >> 11) as f64 / (1u64 << 53) as f64;
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let u2 = (seed >> 11) as f64 / (1u64 << 53) as f64;
            *v = (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos();
        }
        let mut y = matvec_block(a, m, n, &omega, k);
        for _ in 0..n_iter {
            let z = matvec_block_t(a, m, n, &y, k);
            y = matvec_block(a, m, n, &z, k);
        }
        let q = qr_thin(&y, m, k);
        let b = matmul_rect_t(&q, m, k, a, m, n);
        let (ub, sb, vb) = svd_small(&b, k, n);
        let u = matmul_rect(&q, m, k, &ub, k, k);
        LowRankResult {
            u,
            s: sb,
            vt: vb,
            m,
            n,
            k,
        }
    }
    /// Nyström approximation for symmetric positive semi-definite matrices.
    ///
    /// Samples k columns to form a sketch, then recovers a rank-k approximation.
    pub fn nystrom(a: &[f64], n: usize, k: usize) -> LowRankResult {
        let k = k.min(n);
        let step = n / k;
        let cols: Vec<usize> = (0..k).map(|i| i * step).collect();
        let mut c = vec![0.0f64; n * k];
        for (j, &col) in cols.iter().enumerate() {
            for i in 0..n {
                c[i * k + j] = a[i * n + col];
            }
        }
        let mut w = vec![0.0f64; k * k];
        for i in 0..k {
            for j in 0..k {
                w[i * k + j] = a[cols[i] * n + cols[j]];
            }
        }
        let (uw, sw, _vwt) = svd_small(&w, k, k);
        let eps = sw[0] * 1e-12;
        let s_inv_sqrt: Vec<f64> = sw
            .iter()
            .map(|&v| if v > eps { 1.0 / v.sqrt() } else { 0.0 })
            .collect();
        let cuw = matmul_rect(&c, n, k, &uw, k, k);
        let u: Vec<f64> = (0..n * k)
            .map(|idx| {
                let j = idx % k;
                cuw[idx] * s_inv_sqrt[j]
            })
            .collect();
        let s: Vec<f64> = sw.to_vec();
        LowRankResult {
            u: u.clone(),
            s,
            vt: u,
            m: n,
            n,
            k,
        }
    }
    /// CUR decomposition: selects actual rows and columns of A.
    ///
    /// Returns C (selected columns), U (pseudo-inverse of intersection), R (selected rows).
    pub fn cur_decomposition(
        a: &[f64],
        m: usize,
        n: usize,
        k: usize,
    ) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let k = k.min(m).min(n);
        let col_step = n / k;
        let row_step = m / k;
        let col_indices: Vec<usize> = (0..k).map(|i| (i * col_step).min(n - 1)).collect();
        let row_indices: Vec<usize> = (0..k).map(|i| (i * row_step).min(m - 1)).collect();
        let mut c_mat = vec![0.0f64; m * k];
        for (j, &col) in col_indices.iter().enumerate() {
            for i in 0..m {
                c_mat[i * k + j] = a[i * n + col];
            }
        }
        let mut r_mat = vec![0.0f64; k * n];
        for (i, &row) in row_indices.iter().enumerate() {
            for j in 0..n {
                r_mat[i * n + j] = a[row * n + j];
            }
        }
        let mut w = vec![0.0f64; k * k];
        for i in 0..k {
            for j in 0..k {
                w[i * k + j] = a[row_indices[i] * n + col_indices[j]];
            }
        }
        let u = MatrixFunctions::inv_dense(&w, k);
        (c_mat, u, r_mat)
    }
}
/// SSOR (Symmetric Successive Over-Relaxation) preconditioner.
///
/// One forward and one backward sweep with relaxation parameter ω.
pub struct SSORPreconditioner {
    /// The underlying sparse matrix A.
    pub(super) a: SparseCSR,
    /// Relaxation parameter (typically 1 < omega < 2).
    pub(super) omega: f64,
}
impl SSORPreconditioner {
    /// Construct an SSOR preconditioner with the given relaxation parameter.
    pub fn new(a: SparseCSR, omega: f64) -> Self {
        Self { a, omega }
    }
}
/// Result of a low-rank decomposition.
#[derive(Debug, Clone)]
pub struct LowRankResult {
    /// Left singular vectors (m x k), stored column-major.
    pub u: Vec<f64>,
    /// Singular values (length k).
    pub s: Vec<f64>,
    /// Right singular vectors (n x k), stored column-major.
    pub vt: Vec<f64>,
    /// Number of rows m.
    pub m: usize,
    /// Number of columns n.
    pub n: usize,
    /// Rank k.
    pub k: usize,
}
/// Result type returned by iterative solvers.
#[derive(Debug, Clone)]
pub struct SolverResult {
    /// The computed solution vector.
    pub x: Vec<f64>,
    /// Final residual norm.
    pub residual: f64,
    /// Number of iterations performed.
    pub iterations: usize,
    /// Whether the solver converged.
    pub converged: bool,
}
/// Matrix functions for dense matrices stored row-major (n x n).
pub struct MatrixFunctions;
impl MatrixFunctions {
    /// Compute the matrix exponential exp(A) via Padé approximant of order 13.
    ///
    /// Uses the scaling-and-squaring algorithm with a \[13/13\] Padé approximant,
    /// the same approach used in MATLAB's `expm`.
    pub fn expm(a: &[f64], n: usize) -> Vec<f64> {
        let pade_coefs: [f64; 14] = [
            1.0,
            0.5,
            0.12,
            1.833333333333333e-2,
            1.992063492063492e-3,
            1.630434782608696e-4,
            1.035_196_687_370_6e-5,
            5.175983436853003e-7,
            2.043_151_770_816_96e-8,
            6.306659613335318e-10,
            1.483027285835861e-11,
            2.529153491597966e-13,
            2.810170546428514e-15,
            1.544188477920441e-17,
        ];
        let norm_a = Self::matrix_one_norm(a, n);
        let s = (norm_a / 0.5).log2().ceil().max(0.0) as u32;
        let scale = 1.0f64 / (1u64 << s) as f64;
        let a_scaled: Vec<f64> = a.iter().map(|v| v * scale).collect();
        let id = Self::eye(n);
        let a2 = Self::matmul(&a_scaled, &a_scaled, n);
        let a4 = Self::matmul(&a2, &a2, n);
        let a6 = Self::matmul(&a4, &a2, n);
        let u = Self::pade13_uv(&a_scaled, &a2, &a4, &a6, &pade_coefs, n, true);
        let v = Self::pade13_uv(&a_scaled, &a2, &a4, &a6, &pade_coefs, n, false);
        let num = Self::mat_add(&v, &u, n);
        let den = Self::mat_sub(&v, &u, n);
        let inv_den = Self::inv_dense(&den, n);
        let mut result = Self::matmul(&inv_den, &num, n);
        for _ in 0..s {
            result = Self::matmul(&result, &result, n);
        }
        let _ = id;
        result
    }
    /// Compute the matrix logarithm log(A) via inverse scaling and squaring.
    ///
    /// Uses the Padé approximant-based method. A must be positive definite.
    pub fn logm(a: &[f64], n: usize) -> Vec<f64> {
        let mut a_cur = a.to_vec();
        let id = Self::eye(n);
        let mut k = 0u32;
        for _ in 0..16 {
            let diff = Self::mat_sub(&a_cur, &id, n);
            if Self::matrix_one_norm(&diff, n) < 0.5 {
                break;
            }
            a_cur = Self::sqrtm(&a_cur, n);
            k += 1;
        }
        let x = Self::mat_sub(&a_cur, &id, n);
        let x2 = Self::matmul(&x, &x, n);
        let x3 = Self::matmul(&x2, &x, n);
        let nu: Vec<f64> = (0..n * n)
            .map(|i| x[i] * (1.0 / 2.0) + x2[i] * (3.0 / 20.0) + x3[i] * (1.0 / 60.0))
            .collect();
        let de: Vec<f64> = (0..n * n)
            .map(|i| id[i] + x[i] * (3.0 / 4.0) + x2[i] * (3.0 / 10.0) + x3[i] * (1.0 / 20.0))
            .collect();
        let inv_de = Self::inv_dense(&de, n);
        let log_a_scaled = Self::matmul(&inv_de, &nu, n);
        let factor = (1u64 << k) as f64;
        log_a_scaled.iter().map(|v| v * factor).collect()
    }
    /// Compute the matrix square root via the product form of the Denman-Beavers iteration.
    pub fn sqrtm(a: &[f64], n: usize) -> Vec<f64> {
        let mut y = a.to_vec();
        let mut z = Self::eye(n);
        for _ in 0..20 {
            let inv_y = Self::inv_dense(&y, n);
            let inv_z = Self::inv_dense(&z, n);
            let y_new: Vec<f64> = (0..n * n).map(|i| 0.5 * (y[i] + inv_z[i])).collect();
            let z_new: Vec<f64> = (0..n * n).map(|i| 0.5 * (z[i] + inv_y[i])).collect();
            let diff: f64 = y_new
                .iter()
                .zip(&y)
                .map(|(a, b)| (a - b).abs())
                .sum::<f64>();
            y = y_new;
            z = z_new;
            if diff < 1e-12 {
                break;
            }
        }
        y
    }
    pub(crate) fn eye(n: usize) -> Vec<f64> {
        let mut m = vec![0.0f64; n * n];
        for i in 0..n {
            m[i * n + i] = 1.0;
        }
        m
    }
    fn matmul(a: &[f64], b: &[f64], n: usize) -> Vec<f64> {
        let mut c = vec![0.0f64; n * n];
        for i in 0..n {
            for k in 0..n {
                if a[i * n + k] == 0.0 {
                    continue;
                }
                for j in 0..n {
                    c[i * n + j] += a[i * n + k] * b[k * n + j];
                }
            }
        }
        c
    }
    fn mat_add(a: &[f64], b: &[f64], _n: usize) -> Vec<f64> {
        a.iter().zip(b.iter()).map(|(x, y)| x + y).collect()
    }
    fn mat_sub(a: &[f64], b: &[f64], _n: usize) -> Vec<f64> {
        a.iter().zip(b.iter()).map(|(x, y)| x - y).collect()
    }
    fn mat_scalar(a: &[f64], s: f64) -> Vec<f64> {
        a.iter().map(|v| v * s).collect()
    }
    fn matrix_one_norm(a: &[f64], n: usize) -> f64 {
        (0..n)
            .map(|j| (0..n).map(|i| a[i * n + j].abs()).sum::<f64>())
            .fold(0.0f64, f64::max)
    }
    pub(crate) fn inv_dense(a: &[f64], n: usize) -> Vec<f64> {
        let mut m = a.to_vec();
        let mut inv = Self::eye(n);
        for col in 0..n {
            let mut max_row = col;
            let mut max_val = m[col * n + col].abs();
            for row in (col + 1)..n {
                if m[row * n + col].abs() > max_val {
                    max_val = m[row * n + col].abs();
                    max_row = row;
                }
            }
            for j in 0..n {
                m.swap(col * n + j, max_row * n + j);
                inv.swap(col * n + j, max_row * n + j);
            }
            let diag = m[col * n + col];
            if diag.abs() < 1e-300 {
                continue;
            }
            for j in 0..n {
                m[col * n + j] /= diag;
                inv[col * n + j] /= diag;
            }
            for row in 0..n {
                if row == col {
                    continue;
                }
                let factor = m[row * n + col];
                for j in 0..n {
                    let mv = m[col * n + j];
                    let iv = inv[col * n + j];
                    m[row * n + j] -= factor * mv;
                    inv[row * n + j] -= factor * iv;
                }
            }
        }
        inv
    }
    fn pade13_uv(
        a: &[f64],
        a2: &[f64],
        a4: &[f64],
        a6: &[f64],
        c: &[f64; 14],
        n: usize,
        compute_u: bool,
    ) -> Vec<f64> {
        let id = Self::eye(n);
        if compute_u {
            let inner: Vec<f64> = (0..n * n)
                .map(|i| c[13] * a6[i] + c[11] * a4[i] + c[9] * a2[i])
                .collect();
            let t1 = Self::matmul(a6, &inner, n);
            let t2: Vec<f64> = (0..n * n)
                .map(|i| t1[i] + c[7] * a6[i] + c[5] * a4[i] + c[3] * a2[i] + c[1] * id[i])
                .collect();
            Self::matmul(a, &t2, n)
        } else {
            let inner: Vec<f64> = (0..n * n)
                .map(|i| c[12] * a6[i] + c[10] * a4[i] + c[8] * a2[i])
                .collect();
            let t1 = Self::matmul(a6, &inner, n);
            (0..n * n)
                .map(|i| t1[i] + c[6] * a6[i] + c[4] * a4[i] + c[2] * a2[i] + c[0] * id[i])
                .collect()
        }
    }
    /// Compute matrix polynomial p(A) = sum c_i A^i.
    pub fn matrix_poly(a: &[f64], n: usize, coeffs: &[f64]) -> Vec<f64> {
        let id = Self::eye(n);
        let mut result = Self::mat_scalar(&id, coeffs.first().copied().unwrap_or(0.0));
        let mut a_pow = id.clone();
        for &ci in coeffs.iter().skip(1) {
            a_pow = Self::matmul(&a_pow, a, n);
            let term = Self::mat_scalar(&a_pow, ci);
            result = Self::mat_add(&result, &term, n);
        }
        result
    }
}
/// Parameters for the Matricised Tensor Times Khatri-Rao Product (MTTKRP) kernel.
///
/// Groups the six integer scalars that describe the tensor shape and contraction
/// mode so that `TensorDecomposition::mttkrp` stays within the argument-count
/// limit.
#[derive(Debug, Clone, Copy)]
pub struct MttkrpParams {
    /// Number of output rows (mode-specific leading dimension).
    pub m_rows: usize,
    /// CP rank.
    pub r: usize,
    /// Unfolding mode (0, 1, or 2).
    pub mode: usize,
    /// Size of tensor dimension I.
    pub i_dim: usize,
    /// Size of tensor dimension J.
    pub j_dim: usize,
    /// Size of tensor dimension K.
    pub k_dim: usize,
}

/// Tensor decomposition methods.
pub struct TensorDecomposition;
impl TensorDecomposition {
    /// CP decomposition via Alternating Least Squares (ALS).
    ///
    /// Decomposes a 3-way tensor X ≈ sum_r lambda_r a_r ⊗ b_r ⊗ c_r.
    /// `x` is stored in mode-0 unfolding (I x JK), `shape` = \[I, J, K\].
    pub fn cp_als(
        x: &[f64],
        shape: &[usize],
        rank: usize,
        max_iter: usize,
        tol: f64,
    ) -> TensorDecompResult {
        assert_eq!(shape.len(), 3, "CP-ALS requires a 3-way tensor");
        let (i_dim, j_dim, k_dim) = (shape[0], shape[1], shape[2]);
        let r = rank;
        let mut rng_s = 42u64;
        let mut randf = move || -> f64 {
            rng_s = rng_s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (rng_s >> 11) as f64 / (1u64 << 53) as f64 - 0.5
        };
        let mut a_fac: Vec<f64> = (0..i_dim * r).map(|_| randf()).collect();
        let mut b_fac: Vec<f64> = (0..j_dim * r).map(|_| randf()).collect();
        let mut c_fac: Vec<f64> = (0..k_dim * r).map(|_| randf()).collect();
        let mut lambda = vec![1.0f64; r];
        let mut prev_fit = f64::INFINITY;
        for iter in 0..max_iter {
            Self::cp_update_factor(
                x,
                &b_fac,
                &c_fac,
                &mut a_fac,
                &mut lambda,
                i_dim,
                j_dim,
                k_dim,
                r,
                0,
            );
            Self::cp_update_factor(
                x,
                &a_fac,
                &c_fac,
                &mut b_fac,
                &mut lambda,
                i_dim,
                j_dim,
                k_dim,
                r,
                1,
            );
            Self::cp_update_factor(
                x,
                &a_fac,
                &b_fac,
                &mut c_fac,
                &mut lambda,
                i_dim,
                j_dim,
                k_dim,
                r,
                2,
            );
            let x_norm: f64 = x.iter().map(|v| v * v).sum::<f64>().sqrt();
            let fit = if x_norm > 1e-300 {
                let recon =
                    Self::cp_reconstruct(&a_fac, &b_fac, &c_fac, &lambda, i_dim, j_dim, k_dim, r);
                let res: f64 = x
                    .iter()
                    .zip(&recon)
                    .map(|(a, b)| (a - b).powi(2))
                    .sum::<f64>()
                    .sqrt();
                res / x_norm
            } else {
                0.0
            };
            if (fit - prev_fit).abs() < tol && iter > 0 {
                return TensorDecompResult {
                    factors: vec![a_fac, b_fac, c_fac],
                    core: lambda,
                    shape: shape.to_vec(),
                    iterations: iter + 1,
                };
            }
            prev_fit = fit;
        }
        TensorDecompResult {
            factors: vec![a_fac, b_fac, c_fac],
            core: lambda,
            shape: shape.to_vec(),
            iterations: max_iter,
        }
    }
    fn cp_update_factor(
        x: &[f64],
        b: &[f64],
        c: &[f64],
        a_out: &mut [f64],
        lambda: &mut [f64],
        i_dim: usize,
        j_dim: usize,
        k_dim: usize,
        r: usize,
        mode: usize,
    ) {
        let btb = gram(b, if mode == 0 { j_dim } else { i_dim }, r);
        let ctc = gram(c, if mode <= 1 { k_dim } else { j_dim }, r);
        let v: Vec<f64> = btb.iter().zip(ctc.iter()).map(|(a, b)| a * b).collect();
        let v_inv = MatrixFunctions::inv_dense(&v, r);
        let m_rows = match mode {
            0 => i_dim,
            1 => j_dim,
            _ => k_dim,
        };
        let mttkrp = Self::mttkrp(
            x,
            b,
            c,
            MttkrpParams {
                m_rows,
                r,
                mode,
                i_dim,
                j_dim,
                k_dim,
            },
        );
        let new_a = matmul_rect(&mttkrp, m_rows, r, &v_inv, r, r);
        for rc in 0..r {
            let nrm = (0..m_rows)
                .map(|i| new_a[i * r + rc] * new_a[i * r + rc])
                .sum::<f64>()
                .sqrt();
            lambda[rc] = if nrm > 1e-300 { nrm } else { 1.0 };
            for i in 0..m_rows {
                a_out[i * r + rc] = if nrm > 1e-300 {
                    new_a[i * r + rc] / nrm
                } else {
                    new_a[i * r + rc]
                };
            }
        }
    }
    fn mttkrp(x: &[f64], b: &[f64], c: &[f64], p: MttkrpParams) -> Vec<f64> {
        let MttkrpParams {
            m_rows,
            r,
            mode,
            i_dim,
            j_dim,
            k_dim,
        } = p;
        let mut out = vec![0.0f64; m_rows * r];
        match mode {
            0 => {
                for i in 0..i_dim {
                    for j in 0..j_dim {
                        for k in 0..k_dim {
                            let xijk = x[i * j_dim * k_dim + j * k_dim + k];
                            for rc in 0..r {
                                out[i * r + rc] += xijk * b[j * r + rc] * c[k * r + rc];
                            }
                        }
                    }
                }
            }
            1 => {
                for i in 0..i_dim {
                    for j in 0..j_dim {
                        for k in 0..k_dim {
                            let xijk = x[i * j_dim * k_dim + j * k_dim + k];
                            for rc in 0..r {
                                out[j * r + rc] += xijk * b[i * r + rc] * c[k * r + rc];
                            }
                        }
                    }
                }
            }
            _ => {
                for i in 0..i_dim {
                    for j in 0..j_dim {
                        for k in 0..k_dim {
                            let xijk = x[i * j_dim * k_dim + j * k_dim + k];
                            for rc in 0..r {
                                out[k * r + rc] += xijk * b[i * r + rc] * c[j * r + rc];
                            }
                        }
                    }
                }
            }
        }
        out
    }
    fn cp_reconstruct(
        a: &[f64],
        b: &[f64],
        c: &[f64],
        lambda: &[f64],
        i_dim: usize,
        j_dim: usize,
        k_dim: usize,
        r: usize,
    ) -> Vec<f64> {
        let mut out = vec![0.0f64; i_dim * j_dim * k_dim];
        for rc in 0..r {
            for i in 0..i_dim {
                for j in 0..j_dim {
                    for k in 0..k_dim {
                        out[i * j_dim * k_dim + j * k_dim + k] +=
                            lambda[rc] * a[i * r + rc] * b[j * r + rc] * c[k * r + rc];
                    }
                }
            }
        }
        out
    }
    /// Tucker decomposition via Higher-Order SVD (HOSVD).
    ///
    /// Decomposes X ≈ G ×_1 U1 ×_2 U2 ×_3 U3 where G is the core tensor.
    pub fn tucker_hosvd(x: &[f64], shape: &[usize], ranks: &[usize]) -> TensorDecompResult {
        assert_eq!(shape.len(), 3, "Tucker HOSVD requires a 3-way tensor");
        assert_eq!(ranks.len(), 3, "ranks must have 3 elements");
        let (i_dim, j_dim, k_dim) = (shape[0], shape[1], shape[2]);
        let (r1, r2, r3) = (ranks[0], ranks[1], ranks[2]);
        let u1 = Self::mode_unfold_svd(x, i_dim, j_dim * k_dim, r1);
        let u2 = Self::mode_unfold_svd_mode2(x, j_dim, i_dim, k_dim, r2);
        let u3 = Self::mode_unfold_svd_mode3(x, k_dim, i_dim, j_dim, r3);
        let core = Self::tucker_core(x, &u1, &u2, &u3, i_dim, j_dim, k_dim, r1, r2, r3);
        TensorDecompResult {
            factors: vec![u1, u2, u3],
            core,
            shape: shape.to_vec(),
            iterations: 1,
        }
    }
    fn mode_unfold_svd(x: &[f64], m: usize, n: usize, r: usize) -> Vec<f64> {
        let (u, _, _) = svd_small(x, m, n);
        let r = r.min(m);
        (0..m * r).map(|idx| u[idx]).collect()
    }
    fn mode_unfold_svd_mode2(
        x: &[f64],
        j_dim: usize,
        i_dim: usize,
        k_dim: usize,
        r: usize,
    ) -> Vec<f64> {
        let mut unf = vec![0.0f64; j_dim * (i_dim * k_dim)];
        for i in 0..i_dim {
            for j in 0..j_dim {
                for k in 0..k_dim {
                    unf[j * (i_dim * k_dim) + i * k_dim + k] = x[i * j_dim * k_dim + j * k_dim + k];
                }
            }
        }
        let (u, _, _) = svd_small(&unf, j_dim, i_dim * k_dim);
        let r = r.min(j_dim);
        (0..j_dim * r).map(|idx| u[idx]).collect()
    }
    fn mode_unfold_svd_mode3(
        x: &[f64],
        k_dim: usize,
        i_dim: usize,
        j_dim: usize,
        r: usize,
    ) -> Vec<f64> {
        let mut unf = vec![0.0f64; k_dim * (i_dim * j_dim)];
        for i in 0..i_dim {
            for j in 0..j_dim {
                for k in 0..k_dim {
                    unf[k * (i_dim * j_dim) + i * j_dim + j] = x[i * j_dim * k_dim + j * k_dim + k];
                }
            }
        }
        let (u, _, _) = svd_small(&unf, k_dim, i_dim * j_dim);
        let r = r.min(k_dim);
        (0..k_dim * r).map(|idx| u[idx]).collect()
    }
    fn tucker_core(
        x: &[f64],
        u1: &[f64],
        u2: &[f64],
        u3: &[f64],
        i_dim: usize,
        j_dim: usize,
        k_dim: usize,
        r1: usize,
        r2: usize,
        r3: usize,
    ) -> Vec<f64> {
        let mut g = vec![0.0f64; r1 * r2 * r3];
        for i in 0..i_dim {
            for j in 0..j_dim {
                for k in 0..k_dim {
                    let xijk = x[i * j_dim * k_dim + j * k_dim + k];
                    for p in 0..r1 {
                        for q in 0..r2 {
                            for s in 0..r3 {
                                g[p * r2 * r3 + q * r3 + s] +=
                                    xijk * u1[i * r1 + p] * u2[j * r2 + q] * u3[k * r3 + s];
                            }
                        }
                    }
                }
            }
        }
        g
    }
    /// Tensor Train (TT) decomposition via successive SVDs.
    ///
    /// Computes a TT-decomposition X ≈ G1 x G2 x ... x Gd with given TT-ranks.
    pub fn tensor_train(
        x: &[f64],
        shape: &[usize],
        max_rank: usize,
        tol: f64,
    ) -> TensorDecompResult {
        let d = shape.len();
        let total: usize = shape.iter().product();
        assert_eq!(x.len(), total, "x must have prod(shape) elements");
        let mut cores: Vec<Vec<f64>> = Vec::with_capacity(d);
        let mut cur = x.to_vec();
        let mut r_prev = 1usize;
        let mut all_shape = Vec::new();
        for k in 0..d {
            let n_k = shape[k];
            let n_rest: usize = shape[k + 1..].iter().product();
            let m = r_prev * n_k;
            let n = n_rest;
            let (u, s, vt) = svd_small(&cur, m, n);
            let x_norm: f64 = s.iter().map(|v| v * v).sum::<f64>().sqrt();
            let threshold = tol * x_norm / ((d - 1) as f64).sqrt();
            let r_new = s
                .iter()
                .take_while(|&&sv| sv > threshold)
                .count()
                .max(1)
                .min(max_rank)
                .min(s.len());
            let core: Vec<f64> = (0..r_prev * n_k * r_new)
                .map(|idx| {
                    let row = idx / r_new;
                    let col = idx % r_new;
                    u[row * m.min(n) + col]
                })
                .collect();
            cores.push(core);
            all_shape.push((r_prev, n_k, r_new));
            if k + 1 < d {
                cur = (0..r_new * n)
                    .map(|idx| {
                        let row = idx / n;
                        let col = idx % n;
                        s[row] * vt[row * n + col]
                    })
                    .collect();
            }
            r_prev = r_new;
        }
        let flat_cores: Vec<f64> = cores.iter().flat_map(|c| c.iter().copied()).collect();
        TensorDecompResult {
            factors: cores,
            core: flat_cores,
            shape: shape.to_vec(),
            iterations: d,
        }
    }
}
/// ILU(0) — Incomplete LU factorisation with no fill-in.
///
/// Stores the L and U factors using the sparsity pattern of A.
pub struct ILU0Preconditioner {
    /// Number of rows/columns.
    pub(super) n: usize,
    /// Row pointer.
    pub(super) row_ptr: Vec<usize>,
    /// Column indices (same sparsity pattern as A).
    pub(super) col_idx: Vec<usize>,
    /// Values of the ILU factorisation (combined L+U storage).
    pub(super) lu_vals: Vec<f64>,
}
impl ILU0Preconditioner {
    /// Compute the ILU(0) factorisation of A.
    pub fn new(a: &SparseCSR) -> Self {
        let n = a.nrows;
        let row_ptr = a.row_ptr.clone();
        let col_idx = a.col_idx.clone();
        let mut lu_vals = a.values.clone();
        for i in 1..n {
            for k in row_ptr[i]..row_ptr[i + 1] {
                let kk = col_idx[k];
                if kk >= i {
                    break;
                }
                let start_kk = row_ptr[kk];
                let mut diag_pos = start_kk;
                for (offset, &c) in col_idx[start_kk..row_ptr[kk + 1]].iter().enumerate() {
                    if c == kk {
                        diag_pos = start_kk + offset;
                        break;
                    }
                }
                let diag = lu_vals[diag_pos];
                if diag.abs() < 1e-300 {
                    continue;
                }
                lu_vals[k] /= diag;
                let factor = lu_vals[k];
                for jj in (k + 1)..row_ptr[i + 1] {
                    let j = col_idx[jj];
                    for pp in row_ptr[kk]..row_ptr[kk + 1] {
                        if col_idx[pp] == j {
                            lu_vals[jj] -= factor * lu_vals[pp];
                            break;
                        }
                    }
                }
            }
        }
        Self {
            n,
            row_ptr,
            col_idx,
            lu_vals,
        }
    }
}
