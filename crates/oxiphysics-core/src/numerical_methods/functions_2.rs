//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::cholesky;
use super::types::{BifurcationPoint, BifurcationType};

/// Pseudo-arclength continuation for parameter-dependent systems.
///
/// Traces solution branches of F(x, λ) = 0.
/// Returns (x, λ) pairs along the branch.
pub fn pseudo_arclength_continuation<F>(
    f: F,
    x0: f64,
    lambda0: f64,
    ds: f64,
    n_steps: usize,
    tol: f64,
) -> Vec<(f64, f64)>
where
    F: Fn(f64, f64) -> f64,
{
    let mut x = x0;
    let mut lam = lambda0;
    let mut path = vec![(x, lam)];
    let mut tx = 0.0_f64;
    let mut tl = 1.0_f64;
    for _ in 0..n_steps {
        let xp = x + ds * tx;
        let lp = lam + ds * tl;
        let mut xc = xp;
        let mut lc = lp;
        for _ in 0..20 {
            let residual = f(xc, lc);
            let arclength = (xc - xp) * tx + (lc - lp) * tl;
            if residual.abs() < tol && arclength.abs() < tol {
                break;
            }
            let eps = 1e-7;
            let dfdx = (f(xc + eps, lc) - f(xc - eps, lc)) / (2.0 * eps);
            let dfdl = (f(xc, lc + eps) - f(xc, lc - eps)) / (2.0 * eps);
            let det = dfdx * tl - dfdl * tx;
            if det.abs() < 1e-12 {
                break;
            }
            let dx_corr = -(residual * tl - arclength * dfdl) / det;
            let dl_corr = -(dfdx * arclength - dx_corr * tx) / tl.max(1e-12);
            xc += dx_corr * 0.5;
            lc += dl_corr * 0.5;
        }
        let eps = 1e-7;
        let dfdx = (f(xc + eps, lc) - f(xc - eps, lc)) / (2.0 * eps);
        let dfdl = (f(xc, lc + eps) - f(xc, lc - eps)) / (2.0 * eps);
        let tnorm = (dfdx.powi(2) + dfdl.powi(2)).sqrt().max(1e-12);
        tx = -dfdl / tnorm;
        tl = dfdx / tnorm;
        if tx * (xc - x) + tl * (lc - lam) < 0.0 {
            tx = -tx;
            tl = -tl;
        }
        x = xc;
        lam = lc;
        path.push((x, lam));
    }
    path
}
/// Detect bifurcation points on a solution branch.
///
/// Looks for sign changes in the Jacobian determinant (saddle-node) and
/// zero eigenvalues.
pub fn detect_bifurcations<F>(f: F, branch: &[(f64, f64)]) -> Vec<BifurcationPoint>
where
    F: Fn(f64, f64) -> f64,
{
    let mut bifs = Vec::new();
    if branch.len() < 2 {
        return bifs;
    }
    let eps = 1e-6;
    let mut prev_jac = {
        let (x, l) = branch[0];
        (f(x + eps, l) - f(x - eps, l)) / (2.0 * eps)
    };
    for pair in branch.windows(2) {
        let (x0, l0) = pair[0];
        let (x, l) = pair[1];
        let jac = (f(x + eps, l) - f(x - eps, l)) / (2.0 * eps);
        if prev_jac * jac < 0.0 {
            bifs.push(BifurcationPoint {
                x: (x0 + x) / 2.0,
                lambda: (l0 + l) / 2.0,
                bif_type: BifurcationType::SaddleNode,
            });
        }
        prev_jac = jac;
    }
    bifs
}
/// Transpose an m×n matrix.
pub fn transpose(a: &[Vec<f64>]) -> Vec<Vec<f64>> {
    if a.is_empty() {
        return Vec::new();
    }
    let m = a.len();
    let n = a[0].len();
    (0..n).map(|j| (0..m).map(|i| a[i][j]).collect()).collect()
}
/// Matrix-vector multiply.
pub fn mat_vec(a: &[Vec<f64>], x: &[f64]) -> Vec<f64> {
    a.iter()
        .map(|row| row.iter().zip(x.iter()).map(|(aij, xj)| aij * xj).sum())
        .collect()
}
/// Matrix-matrix multiply.
pub fn mat_mul(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let m = a.len();
    if m == 0 || b.is_empty() {
        return Vec::new();
    }
    let n = b[0].len();
    let k = b.len();
    let mut c = vec![vec![0.0f64; n]; m];
    for (i, ci) in c.iter_mut().enumerate() {
        for (j, cij) in ci.iter_mut().enumerate() {
            for l in 0..k {
                *cij += a[i][l] * b[l][j];
            }
        }
    }
    c
}
/// Frobenius norm of a matrix.
pub fn frobenius_norm(a: &[Vec<f64>]) -> f64 {
    a.iter()
        .flat_map(|r| r.iter())
        .map(|x| x.powi(2))
        .sum::<f64>()
        .sqrt()
}
/// Infinity norm (max row sum).
pub fn infinity_norm(a: &[Vec<f64>]) -> f64 {
    a.iter()
        .map(|row| row.iter().map(|x| x.abs()).sum::<f64>())
        .fold(0.0_f64, f64::max)
}
/// Build identity matrix.
pub fn eye(n: usize) -> Vec<Vec<f64>> {
    (0..n)
        .map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
        .collect()
}
/// Reduces a real symmetric matrix to symmetric tridiagonal form T = Qᵀ A Q
/// using Householder reflections (tred2 algorithm, Numerical Recipes §11.3).
///
/// Returns `(d, e, q)` where:
/// - `d` is the diagonal of T (length n)
/// - `e` is the subdiagonal of T (length n, `e[0]` unused)
/// - `q` is the accumulated orthogonal matrix (n×n) stored as q[row][col]
fn householder_tridiag(a: &[Vec<f64>]) -> (Vec<f64>, Vec<f64>, Vec<Vec<f64>>) {
    let n = a.len();
    let mut d = vec![0.0_f64; n];
    let mut e = vec![0.0_f64; n];
    let mut q: Vec<Vec<f64>> = a.iter().map(|row| row.to_vec()).collect();
    for i in (1..n).rev() {
        let l = i;
        let mut h = 0.0_f64;
        let mut scale = 0.0_f64;
        if l > 1 {
            for &qi_k in q[i][..l].iter() {
                scale += qi_k.abs();
            }
            if scale == 0.0 {
                e[i] = q[i][l - 1];
            } else {
                for qi_k in q[i][..l].iter_mut() {
                    *qi_k /= scale;
                    h += *qi_k * *qi_k;
                }
                let f = q[i][l - 1];
                let g = if f >= 0.0 { -h.sqrt() } else { h.sqrt() };
                e[i] = scale * g;
                h -= f * g;
                q[i][l - 1] = f - g;
                let mut f2 = 0.0_f64;
                for j in 0..l {
                    q[j][i] = q[i][j] / h;
                    let mut g2 = 0.0_f64;
                    for (qjk, qik) in q[j][..=j].iter().zip(q[i][..=j].iter()) {
                        g2 += qjk * qik;
                    }
                    for (offset, row_k) in q[j + 1..l].iter().enumerate() {
                        let k = j + 1 + offset;
                        g2 += row_k[j] * q[i][k];
                    }
                    e[j] = g2 / h;
                    f2 += e[j] * q[i][j];
                }
                let hh = f2 / (h + h);
                for j in 0..l {
                    let fj = q[i][j];
                    let gj = e[j] - hh * fj;
                    e[j] = gj;
                    for k in 0..=j {
                        q[j][k] -= fj * e[k] + gj * q[i][k];
                    }
                }
            }
        } else {
            e[i] = q[i][l - 1];
        }
        d[i] = h;
    }
    d[0] = 0.0;
    e[0] = 0.0;
    for i in 0..n {
        let l = i;
        if d[i] != 0.0 {
            for j in 0..l {
                let mut g = 0.0_f64;
                for (qi_k, row_k) in q[i][..l].iter().zip(q[..l].iter()) {
                    g += qi_k * row_k[j];
                }
                for row_k in q[..l].iter_mut() {
                    row_k[j] -= g * row_k[i];
                }
            }
        }
        d[i] = q[i][i];
        q[i][i] = 1.0;
        q[i][..l].fill(0.0);
        for row_j in q[..l].iter_mut() {
            row_j[i] = 0.0;
        }
    }
    (d, e, q)
}
/// Diagonalizes a symmetric tridiagonal matrix in-place using the implicit QL
/// algorithm with Wilkinson shift (tql2 algorithm, Numerical Recipes §11.3).
///
/// On entry `d` contains the diagonal and `e` the subdiagonal (e[0] unused).
/// On exit `d` holds eigenvalues (sorted ascending) and `q` the eigenvectors.
/// Maximum 30 iterations per eigenvalue; continues best-effort if not converged.
fn implicit_ql(d: &mut [f64], e: &mut [f64], q: &mut [Vec<f64>]) {
    let n = d.len();
    for i in 1..n {
        e[i - 1] = e[i];
    }
    e[n - 1] = 0.0;
    for l in 0..n {
        let mut iter = 0_usize;
        loop {
            let mut m = l;
            while m < n - 1 {
                let dd = d[m].abs() + d[m + 1].abs();
                if (e[m].abs() + dd) == dd {
                    break;
                }
                m += 1;
            }
            if m == l {
                break;
            }
            if iter >= 30 {
                break;
            }
            iter += 1;
            let mut g = (d[l + 1] - d[l]) / (2.0 * e[l]);
            let mut r = (g * g + 1.0).sqrt();
            g = d[m] - d[l] + e[l] / (g + if g >= 0.0 { r.abs() } else { -r.abs() });
            let (mut s, mut c, mut p) = (1.0_f64, 1.0_f64, 0.0_f64);
            let mut i = (m as isize) - 1;
            while i >= l as isize {
                let ii = i as usize;
                let f = s * e[ii];
                let b = c * e[ii];
                r = (f * f + g * g).sqrt();
                e[ii + 1] = r;
                if r == 0.0 {
                    d[ii + 1] -= p;
                    e[m] = 0.0;
                    break;
                }
                s = f / r;
                c = g / r;
                g = d[ii + 1] - p;
                r = (d[ii] - g) * s + 2.0 * c * b;
                p = s * r;
                d[ii + 1] = g + p;
                g = c * r - b;
                for row_k in q.iter_mut() {
                    let fk = row_k[ii + 1];
                    row_k[ii + 1] = s * row_k[ii] + c * fk;
                    row_k[ii] = c * row_k[ii] - s * fk;
                }
                i -= 1;
            }
            if r.abs() < 1e-15 && i >= l as isize {
                continue;
            }
            d[l] -= p;
            e[l] = g;
            e[m] = 0.0;
        }
    }
    for i in 0..n - 1 {
        let mut k = i;
        let mut p = d[i];
        for (j, &dj) in d[i + 1..]
            .iter()
            .enumerate()
            .map(|(idx, v)| (idx + i + 1, v))
        {
            if dj < p {
                k = j;
                p = dj;
            }
        }
        if k != i {
            d.swap(k, i);
            for row in q.iter_mut() {
                row.swap(k, i);
            }
        }
    }
}
/// Compute all eigenvalues and eigenvectors of a real symmetric n×n matrix.
///
/// Returns `(eigenvalues, eigenvectors)` where eigenvalues are sorted ascending
/// and `eigenvectors[i]` is the eigenvector corresponding to `eigenvalues[i]`
/// (column-major: `eigvecs[i][row_index]`).
///
/// Algorithm: Householder tridiagonalisation followed by implicit-shift QL.
pub fn symmetric_eigen_n(a: &[Vec<f64>]) -> (Vec<f64>, Vec<Vec<f64>>) {
    let (mut d, mut e, mut q) = householder_tridiag(a);
    implicit_ql(&mut d, &mut e, &mut q);
    let n = d.len();
    let eigvecs: Vec<Vec<f64>> = (0..n).map(|j| (0..n).map(|i| q[i][j]).collect()).collect();
    (d, eigvecs)
}
/// Invert a lower-triangular matrix by forward substitution.
///
/// Returns L⁻¹ such that L * L⁻¹ = I.
fn lower_triangular_inv(l: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = l.len();
    let mut inv = vec![vec![0.0_f64; n]; n];
    for j in 0..n {
        inv[j][j] = if l[j][j].abs() > 1e-15 {
            1.0 / l[j][j]
        } else {
            0.0
        };
        for i in (j + 1)..n {
            let sum: f64 = (j..i).map(|k| l[i][k] * inv[k][j]).sum();
            inv[i][j] = if l[i][i].abs() > 1e-15 {
                -sum / l[i][i]
            } else {
                0.0
            };
        }
    }
    inv
}
/// Solve the generalized symmetric eigenvalue problem Ax = λBx, where B is
/// symmetric positive definite (e.g. the overlap matrix S in quantum chemistry
/// FC = SCε).
///
/// Uses Löwdin canonical orthogonalization: transforms to standard form
/// F′ = L⁻¹ F L⁻ᵀ where L comes from Cholesky(B), solves the standard
/// problem, then back-transforms eigenvectors to the original basis.
///
/// Returns `None` if B is not positive definite (Cholesky fails).
/// Returns `Some((eigenvalues, eigenvectors))` sorted by ascending eigenvalue.
/// Each `eigenvectors[i]` is the i-th eigenvector in the ORIGINAL basis.
pub fn generalized_symmetric_eigen_n(
    a: &[Vec<f64>],
    b: &[Vec<f64>],
) -> Option<(Vec<f64>, Vec<Vec<f64>>)> {
    let n = a.len();
    let l = cholesky(b)?;
    let l_inv = lower_triangular_inv(&l);
    let mut tmp = vec![vec![0.0_f64; n]; n];
    for (i, (tmpi, l_inv_i)) in tmp.iter_mut().zip(l_inv.iter()).enumerate() {
        for (j, tmpij) in tmpi.iter_mut().enumerate() {
            *tmpij = l_inv_i[..=i]
                .iter()
                .zip(a.iter())
                .map(|(&lik, ak)| lik * ak[j])
                .sum();
        }
    }
    let mut f_prime = vec![vec![0.0_f64; n]; n];
    for (i, fpi) in f_prime.iter_mut().enumerate() {
        for (j, fpij) in fpi.iter_mut().enumerate() {
            *fpij = tmp[i][..=j]
                .iter()
                .zip(l_inv[j][..=j].iter())
                .map(|(&tik, &ljk)| tik * ljk)
                .sum();
        }
    }
    let (evals, c_prime) = symmetric_eigen_n(&f_prime);
    let evecs: Vec<Vec<f64>> = (0..n)
        .map(|j| {
            (0..n)
                .map(|i| (i..n).map(|k| l_inv[k][i] * c_prime[j][k]).sum::<f64>())
                .collect()
        })
        .collect();
    Some((evals, evecs))
}
