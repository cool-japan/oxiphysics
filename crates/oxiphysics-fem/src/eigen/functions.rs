//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
use super::types::ModalAnalysis;

/// Find the dominant eigenvalue and eigenvector of a dense n x n matrix
/// using the power iteration method.
///
/// Returns `Some((eigenvalue, eigenvector))` if converged within `max_iter`
/// iterations, or `None` if the iteration failed to converge or the matrix
/// is degenerate.
#[allow(dead_code)]
pub fn power_iteration(
    a: &[Vec<f64>],
    n: usize,
    max_iter: usize,
    tol: f64,
) -> Option<(f64, Vec<f64>)> {
    if n == 0 {
        return None;
    }
    let mut v: Vec<f64> = (0..n).map(|i| 1.0 + 0.01 * (i as f64)).collect();
    normalize(&mut v);
    let mut eigenval = 0.0_f64;
    for iter in 0..max_iter {
        let w = matvec_dense(a, &v, n);
        let vw: f64 = v.iter().zip(w.iter()).map(|(vi, wi)| vi * wi).sum();
        let vv: f64 = v.iter().map(|vi| vi * vi).sum();
        let eigenval_new = if vv.abs() < 1e-60 { 0.0 } else { vw / vv };
        let mut v_new = w;
        let norm = vec_norm(&v_new);
        if norm < 1e-60 {
            return None;
        }
        for x in v_new.iter_mut() {
            *x /= norm;
        }
        let dot: f64 = v.iter().zip(v_new.iter()).map(|(a, b)| a * b).sum();
        if dot < 0.0 {
            for x in v_new.iter_mut() {
                *x = -*x;
            }
        }
        let change = (eigenval_new - eigenval).abs();
        let rel = change / (eigenval_new.abs() + 1e-60);
        v = v_new;
        eigenval = eigenval_new;
        if rel < tol && iter > 0 {
            return Some((eigenval, v));
        }
    }
    Some((eigenval, v))
}
/// Hotelling deflation: remove contribution of a known eigenpair (lambda, phi)
/// from matrix A.
///
/// After deflation: A <- A - lambda * phi * phi^T
/// This allows the next call to `power_iteration` to find the next
/// dominant eigenvalue.
#[allow(dead_code)]
pub fn deflate(a: &mut Vec<Vec<f64>>, eigenval: f64, eigenvec: &[f64]) {
    let n = eigenvec.len();
    for i in 0..n {
        for j in 0..n {
            a[i][j] -= eigenval * eigenvec[i] * eigenvec[j];
        }
    }
}
/// Full Jacobi eigendecomposition for a dense symmetric n x n matrix.
///
/// Returns `(eigenvalues, eigenvectors)` where `eigenvectors[i]` is the
/// i-th eigenvector (column), sorted in ascending order of eigenvalue.
///
/// Uses cyclic Jacobi sweeps until off-diagonal elements are below tolerance.
#[allow(dead_code)]
pub fn jacobi_eigen_dense(a: &[Vec<f64>], n: usize) -> (Vec<f64>, Vec<Vec<f64>>) {
    pub(super) const MAX_SWEEPS: usize = 200;
    pub(super) const EPS: f64 = 1e-14;
    let mut mat: Vec<Vec<f64>> = a.to_vec();
    let mut v: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            let mut row = vec![0.0_f64; n];
            row[i] = 1.0;
            row
        })
        .collect();
    for _ in 0..MAX_SWEEPS {
        let mut off = 0.0_f64;
        for i in 0..n {
            for j in 0..n {
                if i != j {
                    off += mat[i][j] * mat[i][j];
                }
            }
        }
        if off < EPS * EPS {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                let apq = mat[p][q];
                if apq.abs() < EPS {
                    continue;
                }
                let app = mat[p][p];
                let aqq = mat[q][q];
                let tau = (aqq - app) / (2.0 * apq);
                let t = if tau >= 0.0 {
                    1.0 / (tau + (1.0 + tau * tau).sqrt())
                } else {
                    1.0 / (tau - (1.0 + tau * tau).sqrt())
                };
                let c = 1.0 / (1.0 + t * t).sqrt();
                let s = t * c;
                mat[p][p] = app - t * apq;
                mat[q][q] = aqq + t * apq;
                mat[p][q] = 0.0;
                mat[q][p] = 0.0;
                for r in 0..n {
                    if r == p || r == q {
                        continue;
                    }
                    let arp = mat[r][p];
                    let arq = mat[r][q];
                    let new_arp = c * arp - s * arq;
                    let new_arq = s * arp + c * arq;
                    mat[r][p] = new_arp;
                    mat[p][r] = new_arp;
                    mat[r][q] = new_arq;
                    mat[q][r] = new_arq;
                }
                for k in 0..n {
                    let vkp = v[k][p];
                    let vkq = v[k][q];
                    v[k][p] = c * vkp - s * vkq;
                    v[k][q] = s * vkp + c * vkq;
                }
            }
        }
    }
    let eigenvalues: Vec<f64> = (0..n).map(|i| mat[i][i]).collect();
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&a, &b| {
        eigenvalues[a]
            .partial_cmp(&eigenvalues[b])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let sorted_eigenvalues: Vec<f64> = idx.iter().map(|&i| eigenvalues[i]).collect();
    let sorted_eigenvectors: Vec<Vec<f64>> = idx
        .iter()
        .map(|&col| (0..n).map(|row| v[row][col]).collect())
        .collect();
    (sorted_eigenvalues, sorted_eigenvectors)
}
/// Solve the generalized eigenvalue problem K*phi = omega^2*M*phi for a lumped
/// diagonal mass matrix.
///
/// Uses the M^{-1/2} * K * M^{-1/2} transformation to reduce to a standard
/// symmetric eigenvalue problem, then applies Jacobi decomposition.
///
/// Returns `num_modes` sorted `(omega^2, mode_shape)` pairs in ascending order
/// of omega^2.
#[allow(dead_code)]
pub fn generalized_eigen_shift(
    k: &[Vec<f64>],
    m_diag: &[f64],
    n: usize,
    num_modes: usize,
) -> Vec<(f64, Vec<f64>)> {
    let m_inv_sqrt: Vec<f64> = m_diag
        .iter()
        .map(|&m| if m > 0.0 { 1.0 / m.sqrt() } else { 0.0 })
        .collect();
    let k_tilde: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            (0..n)
                .map(|j| m_inv_sqrt[i] * k[i][j] * m_inv_sqrt[j])
                .collect()
        })
        .collect();
    let (eigenvalues, eigenvectors) = jacobi_eigen_dense(&k_tilde, n);
    let num = num_modes.min(n);
    eigenvalues
        .into_iter()
        .zip(eigenvectors)
        .take(num)
        .map(|(omega_sq, psi)| {
            let phi: Vec<f64> = psi
                .iter()
                .zip(m_inv_sqrt.iter())
                .map(|(p, m)| p * m)
                .collect();
            let norm = vec_norm(&phi);
            let phi_norm = if norm > 1e-60 {
                phi.iter().map(|x| x / norm).collect()
            } else {
                phi
            };
            (omega_sq, phi_norm)
        })
        .collect()
}
/// Perform modal analysis: solve K*phi = omega^2*M*phi and return natural frequencies
/// and mode shapes.
#[allow(dead_code)]
pub fn modal_analysis(k: &[Vec<f64>], m_diag: &[f64], num_modes: usize) -> ModalAnalysis {
    let n = k.len();
    let pairs = generalized_eigen_shift(k, m_diag, n, num_modes);
    let mut omega_sq: Vec<f64> = Vec::with_capacity(pairs.len());
    let mut frequencies_hz: Vec<f64> = Vec::with_capacity(pairs.len());
    let mut mode_shapes: Vec<Vec<f64>> = Vec::with_capacity(pairs.len());
    for (osq, shape) in pairs {
        let omega = osq.max(0.0).sqrt();
        let freq_hz = omega / (2.0 * std::f64::consts::PI);
        omega_sq.push(osq);
        frequencies_hz.push(freq_hz);
        mode_shapes.push(shape);
    }
    ModalAnalysis {
        frequencies_hz,
        mode_shapes,
        omega_sq,
    }
}
/// Consistent mass matrix for a 2-node bar element with 2 DOFs per node (4x4).
///
/// Formula: `rho * A * L / 6 * [[2,0,1,0\],[0,2,0,1],[1,0,2,0],[0,1,0,2]]`
#[allow(dead_code)]
pub fn consistent_mass_matrix_bar(length: f64, rho: f64, area: f64) -> [[f64; 4]; 4] {
    let scale = rho * area * length / 6.0;
    [
        [2.0 * scale, 0.0, 1.0 * scale, 0.0],
        [0.0, 2.0 * scale, 0.0, 1.0 * scale],
        [1.0 * scale, 0.0, 2.0 * scale, 0.0],
        [0.0, 1.0 * scale, 0.0, 2.0 * scale],
    ]
}
/// Subspace iteration for finding the `p` smallest eigenpairs of a
/// dense symmetric matrix.
///
/// Starts with `p` random initial vectors and iteratively refines them.
/// Returns up to `p` converged `(eigenvalue, eigenvector)` pairs sorted
/// in ascending order.
#[allow(dead_code)]
pub fn subspace_iteration(
    a: &[Vec<f64>],
    n: usize,
    num_modes: usize,
    max_iter: usize,
    tol: f64,
) -> Vec<(f64, Vec<f64>)> {
    let p = num_modes.min(n);
    if p == 0 || n == 0 {
        return Vec::new();
    }
    let mut q: Vec<Vec<f64>> = (0..p)
        .map(|k| {
            let mut v = vec![0.0; n];
            for i in 0..n {
                let seed = ((i + 1) * (k + 1) * 73 + 17) % 1000;
                v[i] = (seed as f64) / 1000.0 - 0.5;
            }
            v
        })
        .collect();
    gram_schmidt(&mut q, n);
    let mut eigenvalues = vec![0.0; p];
    for _iter in 0..max_iter {
        let mut y: Vec<Vec<f64>> = Vec::with_capacity(p);
        for k in 0..p {
            y.push(matvec_dense(a, &q[k], n));
        }
        let mut b_reduced: Vec<Vec<f64>> = vec![vec![0.0; p]; p];
        for i in 0..p {
            for j in 0..p {
                let dot: f64 = q[i].iter().zip(y[j].iter()).map(|(a, b)| a * b).sum();
                b_reduced[i][j] = dot;
            }
        }
        let (evals, evecs) = jacobi_eigen_dense(&b_reduced, p);
        let mut converged = true;
        for k in 0..p {
            let change = (evals[k] - eigenvalues[k]).abs();
            let rel = change / (evals[k].abs() + 1e-60);
            if rel > tol {
                converged = false;
            }
        }
        eigenvalues = evals.clone();
        let mut q_new: Vec<Vec<f64>> = vec![vec![0.0; n]; p];
        for k in 0..p {
            for i in 0..n {
                for j in 0..p {
                    q_new[k][i] += y[j][i] * evecs[k][j];
                }
            }
        }
        for k in 0..p {
            let norm = vec_norm(&q_new[k]);
            if norm > 1e-60 {
                for x in q_new[k].iter_mut() {
                    *x /= norm;
                }
            }
        }
        q = q_new;
        if converged && _iter > 0 {
            break;
        }
    }
    eigenvalues.into_iter().zip(q).collect()
}
/// Lanczos algorithm for finding eigenvalues of a symmetric matrix.
///
/// Builds a tridiagonal matrix T from Lanczos vectors, then solves
/// the tridiagonal eigenproblem. Returns up to `num_modes` smallest
/// eigenpairs.
///
/// This implementation includes full reorthogonalization to maintain
/// numerical stability.
#[allow(dead_code)]
pub fn lanczos(
    a: &[Vec<f64>],
    n: usize,
    num_modes: usize,
    max_lanczos_steps: usize,
) -> Vec<(f64, Vec<f64>)> {
    if n == 0 || num_modes == 0 {
        return Vec::new();
    }
    let m = max_lanczos_steps.min(n);
    let mut v: Vec<f64> = (0..n).map(|i| 1.0 + 0.01 * (i as f64)).collect();
    normalize(&mut v);
    let mut lanczos_vecs: Vec<Vec<f64>> = Vec::with_capacity(m);
    lanczos_vecs.push(v.clone());
    let mut alpha_vec: Vec<f64> = Vec::with_capacity(m);
    let mut beta_vec: Vec<f64> = Vec::with_capacity(m);
    let mut v_prev = vec![0.0; n];
    let mut beta_prev = 0.0_f64;
    for j in 0..m {
        let mut w = matvec_dense(a, &lanczos_vecs[j], n);
        let alpha: f64 = lanczos_vecs[j]
            .iter()
            .zip(w.iter())
            .map(|(a, b)| a * b)
            .sum();
        alpha_vec.push(alpha);
        for i in 0..n {
            w[i] -= alpha * lanczos_vecs[j][i] + beta_prev * v_prev[i];
        }
        for k in 0..=j {
            let dot: f64 = w
                .iter()
                .zip(lanczos_vecs[k].iter())
                .map(|(a, b)| a * b)
                .sum();
            for i in 0..n {
                w[i] -= dot * lanczos_vecs[k][i];
            }
        }
        let beta = vec_norm(&w);
        beta_vec.push(beta);
        if beta < 1e-14 || j + 1 >= m {
            break;
        }
        for x in w.iter_mut() {
            *x /= beta;
        }
        v_prev = lanczos_vecs[j].clone();
        beta_prev = beta;
        lanczos_vecs.push(w);
    }
    let k = alpha_vec.len();
    let mut t_mat: Vec<Vec<f64>> = vec![vec![0.0; k]; k];
    for i in 0..k {
        t_mat[i][i] = alpha_vec[i];
        if i + 1 < k && i < beta_vec.len() {
            t_mat[i][i + 1] = beta_vec[i];
            t_mat[i + 1][i] = beta_vec[i];
        }
    }
    let (t_evals, t_evecs) = jacobi_eigen_dense(&t_mat, k);
    let num = num_modes.min(k);
    let mut result = Vec::with_capacity(num);
    for mode in 0..num {
        let mut v_mode = vec![0.0; n];
        let num_vecs = lanczos_vecs.len().min(k);
        for j in 0..num_vecs {
            let coeff = t_evecs[mode][j];
            for i in 0..n {
                v_mode[i] += coeff * lanczos_vecs[j][i];
            }
        }
        let norm = vec_norm(&v_mode);
        if norm > 1e-60 {
            for x in v_mode.iter_mut() {
                *x /= norm;
            }
        }
        result.push((t_evals[mode], v_mode));
    }
    result
}
/// Solve the shifted eigenvalue problem `(A - sigma*I)^{-1} * x = mu * x`
/// where eigenvalues of A near sigma correspond to large `mu = 1/(lambda - sigma)`.
///
/// This is useful for finding interior eigenvalues near a specified shift.
///
/// Returns eigenpairs sorted by proximity to `sigma`.
#[allow(dead_code)]
pub fn shift_invert_power(
    a: &[Vec<f64>],
    n: usize,
    sigma: f64,
    num_modes: usize,
    max_iter: usize,
    tol: f64,
) -> Vec<(f64, Vec<f64>)> {
    if n == 0 || num_modes == 0 {
        return Vec::new();
    }
    let mut a_shifted: Vec<Vec<f64>> = a.to_vec();
    for i in 0..n {
        a_shifted[i][i] -= sigma;
    }
    let lu = lu_decompose(&a_shifted, n);
    let mut results: Vec<(f64, Vec<f64>)> = Vec::new();
    let mut a_work = a.to_vec();
    for _mode in 0..num_modes {
        let mut v: Vec<f64> = (0..n)
            .map(|i| 1.0 + 0.02 * (i as f64) + 0.01 * (_mode as f64))
            .collect();
        normalize(&mut v);
        let mut mu = 0.0_f64;
        for iter in 0..max_iter {
            let w = lu_solve(&lu, &v, n);
            let vw: f64 = v.iter().zip(w.iter()).map(|(a, b)| a * b).sum();
            let ww: f64 = w.iter().map(|x| x * x).sum();
            let mu_new = if vw.abs() < 1e-60 { 0.0 } else { ww / vw };
            let mut v_new = w;
            let norm = vec_norm(&v_new);
            if norm < 1e-60 {
                break;
            }
            for x in v_new.iter_mut() {
                *x /= norm;
            }
            let dot: f64 = v.iter().zip(v_new.iter()).map(|(a, b)| a * b).sum();
            if dot < 0.0 {
                for x in v_new.iter_mut() {
                    *x = -*x;
                }
            }
            let rel_change = (mu_new - mu).abs() / (mu_new.abs() + 1e-60);
            v = v_new;
            mu = mu_new;
            if rel_change < tol && iter > 0 {
                break;
            }
        }
        let lambda = if mu.abs() > 1e-60 {
            sigma + 1.0 / mu
        } else {
            sigma
        };
        results.push((lambda, v.clone()));
        deflate(&mut a_work, lambda, &v);
    }
    results.sort_by(|a, b| {
        let da = (a.0 - sigma).abs();
        let db = (b.0 - sigma).abs();
        da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
    });
    results
}
/// Sort eigenvalues by magnitude (ascending).
#[allow(dead_code)]
pub fn sort_eigenvalues_ascending(
    eigenvalues: &[f64],
    eigenvectors: &[Vec<f64>],
) -> (Vec<f64>, Vec<Vec<f64>>) {
    let mut idx: Vec<usize> = (0..eigenvalues.len()).collect();
    idx.sort_by(|&a, &b| {
        eigenvalues[a]
            .partial_cmp(&eigenvalues[b])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let sorted_vals: Vec<f64> = idx.iter().map(|&i| eigenvalues[i]).collect();
    let sorted_vecs: Vec<Vec<f64>> = idx.iter().map(|&i| eigenvectors[i].clone()).collect();
    (sorted_vals, sorted_vecs)
}
/// Sort eigenvalues by absolute magnitude (ascending).
#[allow(dead_code)]
pub fn sort_eigenvalues_by_magnitude(
    eigenvalues: &[f64],
    eigenvectors: &[Vec<f64>],
) -> (Vec<f64>, Vec<Vec<f64>>) {
    let mut idx: Vec<usize> = (0..eigenvalues.len()).collect();
    idx.sort_by(|&a, &b| {
        eigenvalues[a]
            .abs()
            .partial_cmp(&eigenvalues[b].abs())
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let sorted_vals: Vec<f64> = idx.iter().map(|&i| eigenvalues[i]).collect();
    let sorted_vecs: Vec<Vec<f64>> = idx.iter().map(|&i| eigenvectors[i].clone()).collect();
    (sorted_vals, sorted_vecs)
}
/// Track modes across load steps using the Modal Assurance Criterion (MAC).
///
/// The MAC between two mode shapes `phi_a` and `phi_b` is:
/// ```text
/// MAC = (phi_a^T * phi_b)^2 / ((phi_a^T * phi_a) * (phi_b^T * phi_b))
/// ```
/// A MAC value near 1 indicates the modes are the same; near 0 they are different.
#[allow(dead_code)]
pub fn mac_value(phi_a: &[f64], phi_b: &[f64]) -> f64 {
    let ab: f64 = phi_a.iter().zip(phi_b.iter()).map(|(a, b)| a * b).sum();
    let aa: f64 = phi_a.iter().map(|a| a * a).sum();
    let bb: f64 = phi_b.iter().map(|b| b * b).sum();
    if aa * bb < 1e-60 {
        return 0.0;
    }
    (ab * ab) / (aa * bb)
}
/// Compute the MAC matrix between two sets of mode shapes.
///
/// Returns an `m x n` matrix where `m` is the number of modes in set A
/// and `n` is the number in set B.
#[allow(dead_code)]
pub fn mac_matrix(modes_a: &[Vec<f64>], modes_b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let m = modes_a.len();
    let n = modes_b.len();
    let mut mat = vec![vec![0.0; n]; m];
    for i in 0..m {
        for j in 0..n {
            mat[i][j] = mac_value(&modes_a[i], &modes_b[j]);
        }
    }
    mat
}
/// Track mode correspondences between two load steps using the MAC matrix.
///
/// Returns a mapping: `result[i]` is the index in `modes_new` that best
/// corresponds to mode `i` in `modes_old`.
#[allow(dead_code)]
pub fn track_modes(modes_old: &[Vec<f64>], modes_new: &[Vec<f64>]) -> Vec<usize> {
    let mac = mac_matrix(modes_old, modes_new);
    let m = modes_old.len();
    let n = modes_new.len();
    let mut mapping = vec![0usize; m];
    let mut used = vec![false; n];
    for i in 0..m {
        let mut best_j = 0;
        let mut best_mac = -1.0_f64;
        for j in 0..n {
            if !used[j] && mac[i][j] > best_mac {
                best_mac = mac[i][j];
                best_j = j;
            }
        }
        mapping[i] = best_j;
        if best_j < n {
            used[best_j] = true;
        }
    }
    mapping
}
/// Compute effective modal mass for each mode.
///
/// The effective modal mass indicates how much each mode participates
/// in the total response for excitation in a given direction.
///
/// `phi` - mode shape, `m_diag` - diagonal mass, `direction` - unit direction vector
///
/// Returns `Gamma^2 * m_modal` where `Gamma = phi^T * M * r` and
/// `m_modal = phi^T * M * phi`.
#[allow(dead_code)]
pub fn effective_modal_mass(phi: &[f64], m_diag: &[f64], direction: &[f64]) -> f64 {
    let n = phi.len();
    let mut gamma = 0.0;
    let mut m_modal = 0.0;
    for i in 0..n {
        gamma += phi[i] * m_diag[i] * direction[i];
        m_modal += phi[i] * m_diag[i] * phi[i];
    }
    if m_modal.abs() < 1e-60 {
        return 0.0;
    }
    (gamma * gamma) / m_modal
}
/// Dense matrix-vector product: y = A * x
pub(super) fn matvec_dense(a: &[Vec<f64>], x: &[f64], n: usize) -> Vec<f64> {
    (0..n)
        .map(|i| a[i].iter().zip(x.iter()).map(|(aij, xj)| aij * xj).sum())
        .collect()
}
/// Euclidean norm of a vector.
pub(super) fn vec_norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}
/// Normalize a vector in-place (L2 norm = 1).
pub(super) fn normalize(v: &mut Vec<f64>) {
    let norm = vec_norm(v);
    if norm > 1e-60 {
        for x in v.iter_mut() {
            *x /= norm;
        }
    }
}
/// Gram-Schmidt orthogonalization of a set of vectors.
#[allow(dead_code)]
pub(super) fn gram_schmidt(vectors: &mut Vec<Vec<f64>>, n: usize) {
    let p = vectors.len();
    for i in 0..p {
        for j in 0..i {
            let dot: f64 = vectors[i]
                .iter()
                .zip(vectors[j].iter())
                .map(|(a, b)| a * b)
                .sum();
            for k in 0..n {
                vectors[i][k] -= dot * vectors[j][k];
            }
        }
        let norm = vec_norm(&vectors[i]);
        if norm > 1e-60 {
            for x in vectors[i].iter_mut() {
                *x /= norm;
            }
        }
    }
}
/// Simple LU decomposition (Doolittle, no pivoting) for dense matrix.
/// Returns (L, U) stored in a single matrix (L below diagonal, U on/above diagonal).
#[allow(dead_code)]
pub(super) fn lu_decompose(a: &[Vec<f64>], n: usize) -> Vec<Vec<f64>> {
    let mut lu: Vec<Vec<f64>> = a.to_vec();
    for k in 0..n {
        let pivot = lu[k][k];
        if pivot.abs() < 1e-30 {
            lu[k][k] += 1e-12;
        }
        for i in (k + 1)..n {
            lu[i][k] /= lu[k][k];
            for j in (k + 1)..n {
                lu[i][j] -= lu[i][k] * lu[k][j];
            }
        }
    }
    lu
}
/// Solve LU * x = b given the combined LU matrix.
#[allow(dead_code)]
pub(super) fn lu_solve(lu: &[Vec<f64>], b: &[f64], n: usize) -> Vec<f64> {
    let mut y = b.to_vec();
    for i in 0..n {
        for j in 0..i {
            y[i] -= lu[i][j] * y[j];
        }
    }
    let mut x = y;
    for i in (0..n).rev() {
        for j in (i + 1)..n {
            x[i] -= lu[i][j] * x[j];
        }
        if lu[i][i].abs() > 1e-60 {
            x[i] /= lu[i][i];
        }
    }
    x
}
#[cfg(test)]
mod tests {
    use super::*;

    /// 2-DOF spring system: K = \[\[2,-1\\],\[-1,1\]], M = I
    #[test]
    fn test_2dof_spring_system() {
        let k = vec![vec![2.0, -1.0], vec![-1.0, 1.0]];
        let m_diag = vec![1.0, 1.0];
        let omega_sq_low = (3.0 - 5.0_f64.sqrt()) / 2.0;
        let omega_sq_high = (3.0 + 5.0_f64.sqrt()) / 2.0;
        let pairs = generalized_eigen_shift(&k, &m_diag, 2, 2);
        assert_eq!(pairs.len(), 2, "should return 2 modes");
        let (osq0, _) = &pairs[0];
        let (osq1, _) = &pairs[1];
        assert!(
            (osq0 - omega_sq_low).abs() < 1e-8,
            "omega^2_0 = {osq0}, expected {omega_sq_low}"
        );
        assert!(
            (osq1 - omega_sq_high).abs() < 1e-8,
            "omega^2_1 = {osq1}, expected {omega_sq_high}"
        );
    }
    /// Single DOF: K = \[\[k0\\]], M = [m0]
    #[test]
    fn test_single_dof() {
        let k0 = 400.0_f64;
        let m0 = 4.0_f64;
        let k = vec![vec![k0]];
        let m_diag = vec![m0];
        let result = modal_analysis(&k, &m_diag, 1);
        let expected_omega_sq = k0 / m0;
        let expected_freq_hz = expected_omega_sq.sqrt() / (2.0 * std::f64::consts::PI);
        assert_eq!(result.omega_sq.len(), 1);
        assert!(
            (result.omega_sq[0] - expected_omega_sq).abs() < 1e-8,
            "omega^2 = {}, expected {expected_omega_sq}",
            result.omega_sq[0]
        );
        assert!(
            (result.frequencies_hz[0] - expected_freq_hz).abs() < 1e-8,
            "f_hz = {}, expected {expected_freq_hz}",
            result.frequencies_hz[0]
        );
    }
    /// Frequencies must be sorted ascending.
    #[test]
    fn test_frequencies_sorted() {
        let k = vec![
            vec![9.0, 0.0, 0.0],
            vec![0.0, 4.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];
        let m_diag = vec![1.0, 1.0, 1.0];
        let result = modal_analysis(&k, &m_diag, 3);
        for i in 0..result.frequencies_hz.len() - 1 {
            assert!(
                result.frequencies_hz[i] <= result.frequencies_hz[i + 1] + 1e-10,
                "frequencies not sorted: {} > {}",
                result.frequencies_hz[i],
                result.frequencies_hz[i + 1]
            );
        }
        assert!(
            (result.omega_sq[0] - 1.0).abs() < 1e-8,
            "omega^2_0 = {}",
            result.omega_sq[0]
        );
        assert!(
            (result.omega_sq[1] - 4.0).abs() < 1e-8,
            "omega^2_1 = {}",
            result.omega_sq[1]
        );
        assert!(
            (result.omega_sq[2] - 9.0).abs() < 1e-8,
            "omega^2_2 = {}",
            result.omega_sq[2]
        );
    }
    /// Test power_iteration on a simple 2x2 diagonal matrix.
    #[test]
    fn test_power_iteration_diagonal() {
        let a = vec![vec![5.0, 0.0], vec![0.0, 2.0]];
        let result = power_iteration(&a, 2, 1000, 1e-10);
        assert!(result.is_some(), "should converge");
        let (eigenval, _eigenvec) = result.unwrap();
        assert!(
            (eigenval - 5.0).abs() < 1e-6,
            "dominant eigenvalue = {eigenval}, expected 5.0"
        );
    }
    /// Test Hotelling deflation.
    #[test]
    fn test_deflation() {
        let mut a = vec![vec![5.0, 0.0], vec![0.0, 2.0]];
        let eigenvec = vec![1.0, 0.0];
        deflate(&mut a, 5.0, &eigenvec);
        assert!((a[0][0]).abs() < 1e-12, "a[0][0] = {}", a[0][0]);
        assert!((a[1][1] - 2.0).abs() < 1e-12, "a[1][1] = {}", a[1][1]);
    }
    /// Test consistent mass matrix for a bar element.
    #[test]
    fn test_consistent_mass_matrix_bar() {
        let length = 1.0;
        let rho = 1.0;
        let area = 1.0;
        let m = consistent_mass_matrix_bar(length, rho, area);
        let s = 1.0 / 6.0;
        assert!((m[0][0] - 2.0 * s).abs() < 1e-12, "m[0][0] = {}", m[0][0]);
        assert!((m[0][2] - 1.0 * s).abs() < 1e-12, "m[0][2] = {}", m[0][2]);
        assert!((m[1][1] - 2.0 * s).abs() < 1e-12, "m[1][1] = {}", m[1][1]);
        assert!((m[1][3] - 1.0 * s).abs() < 1e-12, "m[1][3] = {}", m[1][3]);
        assert!((m[0][1]).abs() < 1e-12, "m[0][1] should be 0");
        assert!((m[2][3]).abs() < 1e-12, "m[2][3] should be 0");
    }
    /// Test jacobi_eigen_dense on a known 2x2 symmetric matrix.
    #[test]
    fn test_jacobi_eigen_dense_2x2() {
        let a = vec![vec![3.0, 1.0], vec![1.0, 3.0]];
        let (evals, evecs) = jacobi_eigen_dense(&a, 2);
        assert!((evals[0] - 2.0).abs() < 1e-10, "eval[0] = {}", evals[0]);
        assert!((evals[1] - 4.0).abs() < 1e-10, "eval[1] = {}", evals[1]);
        for i in 0..2 {
            let v = &evecs[i];
            let av = matvec_dense(&a, v, 2);
            for j in 0..2 {
                let diff = (av[j] - evals[i] * v[j]).abs();
                assert!(diff < 1e-10, "residual at ({i},{j}) = {diff}");
            }
        }
    }
    /// Subspace iteration should find eigenvalues of a diagonal matrix.
    #[test]
    fn test_subspace_iteration_diagonal() {
        let a = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 4.0, 0.0],
            vec![0.0, 0.0, 9.0],
        ];
        let pairs = subspace_iteration(&a, 3, 3, 200, 1e-8);
        assert_eq!(pairs.len(), 3);
        let mut evals: Vec<f64> = pairs.iter().map(|(e, _)| *e).collect();
        evals.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((evals[0] - 1.0).abs() < 0.1, "eval[0] = {}", evals[0]);
        assert!((evals[1] - 4.0).abs() < 0.1, "eval[1] = {}", evals[1]);
        assert!((evals[2] - 9.0).abs() < 0.1, "eval[2] = {}", evals[2]);
    }
    /// Lanczos should find eigenvalues of a simple symmetric matrix.
    #[test]
    fn test_lanczos_diagonal() {
        let a = vec![
            vec![2.0, 0.0, 0.0],
            vec![0.0, 5.0, 0.0],
            vec![0.0, 0.0, 8.0],
        ];
        let pairs = lanczos(&a, 3, 3, 3);
        assert!(!pairs.is_empty(), "Lanczos should return results");
        let mut evals: Vec<f64> = pairs.iter().map(|(e, _)| *e).collect();
        evals.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!(
            (evals[0] - 2.0).abs() < 0.5,
            "smallest eigenvalue = {}, expected ~2.0",
            evals[0]
        );
    }
    /// Lanczos on a 2x2 matrix should find both eigenvalues.
    #[test]
    fn test_lanczos_2x2() {
        let a = vec![vec![3.0, 1.0], vec![1.0, 3.0]];
        let pairs = lanczos(&a, 2, 2, 2);
        assert_eq!(pairs.len(), 2);
        let mut evals: Vec<f64> = pairs.iter().map(|(e, _)| *e).collect();
        evals.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((evals[0] - 2.0).abs() < 0.1, "eval[0] = {}", evals[0]);
        assert!((evals[1] - 4.0).abs() < 0.1, "eval[1] = {}", evals[1]);
    }
    /// Shift-invert should find eigenvalues near the shift.
    #[test]
    fn test_shift_invert_near_target() {
        let a = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 5.0, 0.0],
            vec![0.0, 0.0, 10.0],
        ];
        let pairs = shift_invert_power(&a, 3, 4.5, 1, 500, 1e-8);
        assert!(!pairs.is_empty());
        let (lambda, _) = &pairs[0];
        assert!(
            (*lambda - 5.0).abs() < 0.5,
            "eigenvalue near 4.5 = {lambda}, expected ~5.0"
        );
    }
    /// Eigenvalue sorting should produce ascending order.
    #[test]
    fn test_sort_ascending() {
        let evals = vec![9.0, 1.0, 4.0];
        let evecs = vec![
            vec![0.0, 0.0, 1.0],
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
        ];
        let (sorted_vals, _sorted_vecs) = sort_eigenvalues_ascending(&evals, &evecs);
        assert!((sorted_vals[0] - 1.0).abs() < 1e-12);
        assert!((sorted_vals[1] - 4.0).abs() < 1e-12);
        assert!((sorted_vals[2] - 9.0).abs() < 1e-12);
    }
    /// Sorting by magnitude.
    #[test]
    fn test_sort_by_magnitude() {
        let evals = vec![-5.0, 1.0, 3.0];
        let evecs = vec![vec![1.0], vec![1.0], vec![1.0]];
        let (sorted_vals, _) = sort_eigenvalues_by_magnitude(&evals, &evecs);
        assert!((sorted_vals[0] - 1.0).abs() < 1e-12);
        assert!((sorted_vals[1] - 3.0).abs() < 1e-12);
        assert!((sorted_vals[2] - (-5.0)).abs() < 1e-12);
    }
    /// MAC of a mode with itself should be 1.
    #[test]
    fn test_mac_self() {
        let phi = vec![1.0, 0.0, 0.0];
        let m = mac_value(&phi, &phi);
        assert!((m - 1.0).abs() < 1e-12, "MAC self = {m}, expected 1.0");
    }
    /// MAC of orthogonal modes should be 0.
    #[test]
    fn test_mac_orthogonal() {
        let phi_a = vec![1.0, 0.0, 0.0];
        let phi_b = vec![0.0, 1.0, 0.0];
        let m = mac_value(&phi_a, &phi_b);
        assert!(m.abs() < 1e-12, "MAC orthogonal = {m}, expected 0.0");
    }
    /// MAC matrix should be identity for identical mode sets.
    #[test]
    fn test_mac_matrix_identity() {
        let modes = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];
        let mat = mac_matrix(&modes, &modes);
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (mat[i][j] - expected).abs() < 1e-12,
                    "MAC[{i}][{j}] = {}, expected {expected}",
                    mat[i][j]
                );
            }
        }
    }
    /// Mode tracking should correctly identify correspondences.
    #[test]
    fn test_track_modes_identity() {
        let modes = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];
        let mapping = track_modes(&modes, &modes);
        assert_eq!(mapping, vec![0, 1, 2]);
    }
    /// Mode tracking with swapped modes.
    #[test]
    fn test_track_modes_swapped() {
        let modes_old = vec![vec![1.0, 0.0, 0.0], vec![0.0, 1.0, 0.0]];
        let modes_new = vec![vec![0.0, 1.0, 0.0], vec![1.0, 0.0, 0.0]];
        let mapping = track_modes(&modes_old, &modes_new);
        assert_eq!(mapping[0], 1, "old mode 0 should map to new mode 1");
        assert_eq!(mapping[1], 0, "old mode 1 should map to new mode 0");
    }
    /// Effective modal mass for a unit mass system.
    #[test]
    fn test_effective_modal_mass() {
        let phi = vec![1.0, 0.0, 0.0];
        let m_diag = vec![1.0, 1.0, 1.0];
        let direction = vec![1.0, 0.0, 0.0];
        let emm = effective_modal_mass(&phi, &m_diag, &direction);
        assert!(
            (emm - 1.0).abs() < 1e-12,
            "effective modal mass = {emm}, expected 1.0"
        );
    }
    /// Effective modal mass for orthogonal direction should be 0.
    #[test]
    fn test_effective_modal_mass_orthogonal() {
        let phi = vec![1.0, 0.0, 0.0];
        let m_diag = vec![1.0, 1.0, 1.0];
        let direction = vec![0.0, 1.0, 0.0];
        let emm = effective_modal_mass(&phi, &m_diag, &direction);
        assert!(
            emm.abs() < 1e-12,
            "should be 0 for orthogonal direction, got {emm}"
        );
    }
    /// Gram-Schmidt should produce orthonormal vectors.
    #[test]
    fn test_gram_schmidt() {
        let mut vecs = vec![
            vec![1.0, 1.0, 0.0],
            vec![1.0, 0.0, 1.0],
            vec![0.0, 1.0, 1.0],
        ];
        gram_schmidt(&mut vecs, 3);
        for i in 0..3 {
            for j in 0..3 {
                let dot: f64 = vecs[i].iter().zip(vecs[j].iter()).map(|(a, b)| a * b).sum();
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (dot - expected).abs() < 1e-10,
                    "dot(v{i}, v{j}) = {dot}, expected {expected}"
                );
            }
        }
    }
    /// LU decomposition and solve should recover the solution.
    #[test]
    fn test_lu_solve() {
        let a = vec![vec![2.0, 1.0], vec![1.0, 3.0]];
        let b = vec![5.0, 7.0];
        let lu = lu_decompose(&a, 2);
        let x = lu_solve(&lu, &b, 2);
        let ax0 = 2.0 * x[0] + 1.0 * x[1];
        let ax1 = 1.0 * x[0] + 3.0 * x[1];
        assert!((ax0 - 5.0).abs() < 1e-10, "Ax[0] = {ax0}, expected 5");
        assert!((ax1 - 7.0).abs() < 1e-10, "Ax[1] = {ax1}, expected 7");
    }
}
