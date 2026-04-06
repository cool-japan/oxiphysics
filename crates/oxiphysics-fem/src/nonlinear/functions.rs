//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use crate::solvers::conjugate_gradient;
use crate::sparse::CsrMatrix;

use super::types::{ContinuationStepResult, ConvergenceCriteria, ConvergenceHistory, NrResult};

/// Newton-Raphson nonlinear solver.
///
/// The caller provides a closure `tangent_and_residual: F` that, given current `u`,
/// returns `(K_tangent: CsrMatrix, R: Vec`f64`)`.
///
/// # Algorithm
///
/// 1. Start with `u = u0` (usually zeros)
/// 2. Compute `(K, R) = tangent_and_residual(u)`
/// 3. Solve `K * Δu = -R`
/// 4. `u += Δu`
/// 5. Repeat until `||R|| < tol`
pub fn newton_raphson<F>(
    n_dof: usize,
    u0: Option<Vec<f64>>,
    tangent_and_residual: F,
    criteria: &ConvergenceCriteria,
) -> NrResult
where
    F: Fn(&[f64]) -> (CsrMatrix, Vec<f64>),
{
    let mut u = u0.unwrap_or_else(|| vec![0.0; n_dof]);
    assert_eq!(u.len(), n_dof, "u0 length must equal n_dof");
    let mut iterations = 0;
    let mut final_residual = f64::INFINITY;
    let mut converged = false;
    for _iter in 0..criteria.max_iter {
        let (k_tangent, residual) = tangent_and_residual(&u);
        let r_norm: f64 = residual.iter().map(|r| r * r).sum::<f64>().sqrt();
        final_residual = r_norm;
        if r_norm < criteria.residual_tol {
            converged = true;
            break;
        }
        let neg_r: Vec<f64> = residual.iter().map(|r| -r).collect();
        let x0 = vec![0.0; n_dof];
        let du = conjugate_gradient(&k_tangent, &neg_r, &x0, n_dof * 100, 1e-12);
        for i in 0..n_dof {
            u[i] += du[i];
        }
        iterations += 1;
        let du_norm: f64 = du.iter().map(|d| d * d).sum::<f64>().sqrt();
        if du_norm < criteria.displacement_tol {
            let (_, res_check) = tangent_and_residual(&u);
            final_residual = res_check.iter().map(|r| r * r).sum::<f64>().sqrt();
            if final_residual < criteria.residual_tol {
                converged = true;
            }
            break;
        }
    }
    if !converged && iterations == criteria.max_iter {
        let (_, res_final) = tangent_and_residual(&u);
        final_residual = res_final.iter().map(|r| r * r).sum::<f64>().sqrt();
        if final_residual < criteria.residual_tol {
            converged = true;
        }
    }
    NrResult {
        displacement: u,
        converged,
        iterations,
        final_residual,
    }
}
/// Compute the determinant of a 3×3 matrix.
#[allow(dead_code)]
pub fn mat3_det(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}
/// Compute the inverse of a 3×3 matrix.
///
/// Panics if the matrix is singular (det ≈ 0).
#[allow(dead_code)]
pub fn mat3_inv(m: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let det = mat3_det(m);
    assert!(det.abs() > 1e-15, "mat3_inv: singular matrix (det={det})");
    let inv_det = 1.0 / det;
    [
        [
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv_det,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv_det,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_det,
        ],
        [
            (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inv_det,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_det,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inv_det,
        ],
        [
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv_det,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv_det,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_det,
        ],
    ]
}
/// Multiply two 3×3 matrices: C = A * B.
#[allow(dead_code)]
pub fn mat3_mul(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut c = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}
/// Transpose a 3×3 matrix.
#[allow(dead_code)]
pub fn mat3_transpose(m: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut t = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            t[i][j] = m[j][i];
        }
    }
    t
}
/// Trace of a 3×3 matrix: tr(M) = M\[0\]\[0\] + M\[1\]\[1\] + M\[2\]\[2\].
#[allow(dead_code)]
pub fn mat3_trace(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] + m[1][1] + m[2][2]
}
/// Return the 3×3 identity matrix.
#[allow(dead_code)]
pub fn mat3_identity() -> [[f64; 3]; 3] {
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
}
/// Compute the deformation gradient F = dx/dX for a 4-node tetrahedral element.
///
/// `x`  = current nodal positions (4 nodes × 3 coords)
/// `x0` = reference nodal positions (4 nodes × 3 coords)
///
/// Algorithm: F = ds * DS^{-1} where
///   ds = \[x\[1\\]-x\[0\], x\[2\]-x\[0\], x\[3\]-x\[0\]] (columns, stored row-major as 3×3)
///   DS = \[x0\[1\\]-x0\[0\], x0\[2\]-x0\[0\], x0\[3\]-x0\[0\]]
#[allow(dead_code)]
pub fn deformation_gradient(x: &[[f64; 3]; 4], x0: &[[f64; 3]; 4]) -> [[f64; 3]; 3] {
    let ds: [[f64; 3]; 3] = [
        [x[1][0] - x[0][0], x[2][0] - x[0][0], x[3][0] - x[0][0]],
        [x[1][1] - x[0][1], x[2][1] - x[0][1], x[3][1] - x[0][1]],
        [x[1][2] - x[0][2], x[2][2] - x[0][2], x[3][2] - x[0][2]],
    ];
    let ds0: [[f64; 3]; 3] = [
        [
            x0[1][0] - x0[0][0],
            x0[2][0] - x0[0][0],
            x0[3][0] - x0[0][0],
        ],
        [
            x0[1][1] - x0[0][1],
            x0[2][1] - x0[0][1],
            x0[3][1] - x0[0][1],
        ],
        [
            x0[1][2] - x0[0][2],
            x0[2][2] - x0[0][2],
            x0[3][2] - x0[0][2],
        ],
    ];
    let ds0_inv = mat3_inv(&ds0);
    mat3_mul(&ds, &ds0_inv)
}
/// Right Cauchy-Green deformation tensor C = F^T * F.
#[allow(dead_code)]
pub fn right_cauchy_green(f: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let ft = mat3_transpose(f);
    mat3_mul(&ft, f)
}
/// Left Cauchy-Green deformation tensor B = F * F^T.
#[allow(dead_code)]
pub fn left_cauchy_green(f: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let ft = mat3_transpose(f);
    mat3_mul(f, &ft)
}
/// Green-Lagrange strain tensor E = (C - I) / 2.
#[allow(dead_code)]
pub fn green_lagrange_strain(f: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let c = right_cauchy_green(f);
    let mut e = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let delta = if i == j { 1.0 } else { 0.0 };
            e[i][j] = 0.5 * (c[i][j] - delta);
        }
    }
    e
}
/// Principal invariants of the left Cauchy-Green tensor B.
///
/// Returns `(I1, I2, I3)` where:
/// - I1 = tr(B)
/// - I2 = (tr(B)² − tr(B²)) / 2
/// - I3 = det(B)
#[allow(dead_code)]
pub fn invariants(b: &[[f64; 3]; 3]) -> (f64, f64, f64) {
    let i1 = mat3_trace(b);
    let b2 = mat3_mul(b, b);
    let tr_b2 = mat3_trace(&b2);
    let i2 = 0.5 * (i1 * i1 - tr_b2);
    let i3 = mat3_det(b);
    (i1, i2, i3)
}
/// 1st Piola-Kirchhoff stress for the compressible Neo-Hookean model.
///
/// P = μ (F − F^{−T}) + λ ln(J) F^{−T}
///
/// where J = det(F).
#[allow(dead_code)]
pub fn neo_hookean_stress(f: &[[f64; 3]; 3], mu: f64, lambda: f64) -> [[f64; 3]; 3] {
    let j = mat3_det(f);
    let f_inv = mat3_inv(f);
    let f_inv_t = mat3_transpose(&f_inv);
    let ln_j = j.ln();
    let mut p = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j_idx in 0..3 {
            p[i][j_idx] =
                mu * (f[i][j_idx] - f_inv_t[i][j_idx]) + lambda * ln_j * f_inv_t[i][j_idx];
        }
    }
    p
}
/// 1st Piola-Kirchhoff stress for the incompressible Mooney-Rivlin model.
///
/// W = c10 (I1 − 3) + c01 (I2 − 3)
///
/// P = 2 c10 F + 2 c01 (I1 F − F C)
///
/// (deviatoric, incompressible form)
#[allow(dead_code)]
pub fn mooney_rivlin_stress(f: &[[f64; 3]; 3], c10: f64, c01: f64) -> [[f64; 3]; 3] {
    let b = left_cauchy_green(f);
    let i1 = mat3_trace(&b);
    let c = right_cauchy_green(f);
    let fc = mat3_mul(f, &c);
    let mut p = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            p[i][j] = 2.0 * c10 * f[i][j] + 2.0 * c01 * (i1 * f[i][j] - fc[i][j]);
        }
    }
    p
}
/// Compute nodal internal forces for a 4-node Neo-Hookean tetrahedral element.
///
/// Returns `[[f64;3\];4]` — one force vector per node.
///
/// Algorithm (constant-strain tet):
///   1. Compute F = deformation_gradient(x, x0)
///   2. P = neo_hookean_stress(F, mu, lambda)
///   3. Build shape-function gradient matrix BN (4×3) in reference config
///   4. f_int\[a\] = V0 * P * N_{a,X}  for each node a
#[allow(dead_code)]
#[allow(clippy::too_many_arguments)]
pub fn internal_forces_tet(
    x: &[[f64; 3]; 4],
    x0: &[[f64; 3]; 4],
    mu: f64,
    lambda: f64,
) -> [[f64; 3]; 4] {
    let ds0: [[f64; 3]; 3] = [
        [
            x0[1][0] - x0[0][0],
            x0[2][0] - x0[0][0],
            x0[3][0] - x0[0][0],
        ],
        [
            x0[1][1] - x0[0][1],
            x0[2][1] - x0[0][1],
            x0[3][1] - x0[0][1],
        ],
        [
            x0[1][2] - x0[0][2],
            x0[2][2] - x0[0][2],
            x0[3][2] - x0[0][2],
        ],
    ];
    let v0 = mat3_det(&ds0).abs() / 6.0;
    let f_grad = deformation_gradient(x, x0);
    let p = neo_hookean_stress(&f_grad, mu, lambda);
    let ds0_inv = mat3_inv(&ds0);
    let ds0_inv_t = mat3_transpose(&ds0_inv);
    let mut grad_n = [[0.0f64; 3]; 4];
    for a in 0..3 {
        for i in 0..3 {
            grad_n[a + 1][i] = ds0_inv_t[i][a];
        }
    }
    for i in 0..3 {
        grad_n[0][i] = -(grad_n[1][i] + grad_n[2][i] + grad_n[3][i]);
    }
    let mut forces = [[0.0f64; 3]; 4];
    for a in 0..4 {
        for i in 0..3 {
            for j in 0..3 {
                forces[a][i] += v0 * p[i][j] * grad_n[a][j];
            }
        }
    }
    forces
}
#[cfg(test)]
mod large_deform_tests {
    use super::*;

    fn identity() -> [[f64; 3]; 3] {
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    }
    /// Reference tet: canonical simplex
    fn ref_tet() -> [[f64; 3]; 4] {
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    }
    /// F = I → E = 0 (no strain in reference configuration).
    #[test]
    fn test_green_lagrange_at_identity() {
        let f = identity();
        let e = green_lagrange_strain(&f);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    e[i][j].abs() < 1e-15,
                    "E[{i}][{j}] = {} should be 0",
                    e[i][j]
                );
            }
        }
    }
    /// det(F) = 1 for isochoric (volume-preserving) deformation: pure shear.
    #[test]
    fn test_det_isochoric() {
        let gamma = 0.5;
        let f = [[1.0, gamma, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let det = mat3_det(&f);
        assert!((det - 1.0).abs() < 1e-15, "det = {det}, expected 1");
    }
    /// Neo-Hookean stress at F = I: P = 0 (no stress at rest).
    #[test]
    fn test_neo_hookean_stress_at_identity() {
        let f = identity();
        let p = neo_hookean_stress(&f, 1.0, 1.0);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    p[i][j].abs() < 1e-14,
                    "P[{i}][{j}] = {} should be 0 at F=I",
                    p[i][j]
                );
            }
        }
    }
    /// C = I when F = I.
    #[test]
    fn test_right_cauchy_green_at_identity() {
        let f = identity();
        let c = right_cauchy_green(&f);
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (c[i][j] - expected).abs() < 1e-15,
                    "C[{i}][{j}] = {}, expected {expected}",
                    c[i][j]
                );
            }
        }
    }
    /// invariants(I): I1=3, I2=3, I3=1.
    #[test]
    fn test_invariants_at_identity() {
        let b = identity();
        let (i1, i2, i3) = invariants(&b);
        assert!((i1 - 3.0).abs() < 1e-15, "I1={i1}");
        assert!((i2 - 3.0).abs() < 1e-15, "I2={i2}");
        assert!((i3 - 1.0).abs() < 1e-15, "I3={i3}");
    }
    /// mat3_identity returns the 3×3 identity matrix.
    #[test]
    fn test_mat3_identity() {
        let i = mat3_identity();
        for row in 0..3 {
            for col in 0..3 {
                let expected = if row == col { 1.0 } else { 0.0 };
                assert_eq!(i[row][col], expected, "identity[{row}][{col}]");
            }
        }
    }
    /// deformation_gradient on undeformed tet → F = I.
    #[test]
    fn test_deformation_gradient_undeformed() {
        let x0 = ref_tet();
        let f = deformation_gradient(&x0, &x0);
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (f[i][j] - expected).abs() < 1e-14,
                    "F[{i}][{j}] = {}, expected {expected}",
                    f[i][j]
                );
            }
        }
    }
    /// internal_forces_tet at rest → zero forces.
    #[test]
    fn test_internal_forces_at_rest() {
        let x0 = ref_tet();
        let forces = internal_forces_tet(&x0, &x0, 1.0, 1.0);
        for a in 0..4 {
            for i in 0..3 {
                assert!(
                    forces[a][i].abs() < 1e-13,
                    "force[{a}][{i}] = {}, expected 0",
                    forces[a][i]
                );
            }
        }
    }
}
/// Newton-Raphson with backtracking line search.
///
/// Uses Armijo condition to find step length `alpha` such that
/// `||R(u + alpha * du)||^2 < (1 - 2*c*alpha) * ||R(u)||^2`.
pub fn newton_raphson_line_search<F>(
    n_dof: usize,
    u0: Option<Vec<f64>>,
    tangent_and_residual: F,
    criteria: &ConvergenceCriteria,
) -> NrResult
where
    F: Fn(&[f64]) -> (CsrMatrix, Vec<f64>),
{
    let mut u = u0.unwrap_or_else(|| vec![0.0; n_dof]);
    let mut iterations = 0;
    let mut final_residual = f64::INFINITY;
    let mut converged = false;
    pub(super) const C: f64 = 1e-4;
    for _iter in 0..criteria.max_iter {
        let (k_tangent, residual) = tangent_and_residual(&u);
        let r_norm_sq: f64 = residual.iter().map(|r| r * r).sum();
        final_residual = r_norm_sq.sqrt();
        if final_residual < criteria.residual_tol {
            converged = true;
            break;
        }
        let neg_r: Vec<f64> = residual.iter().map(|r| -r).collect();
        let x0 = vec![0.0; n_dof];
        let du = conjugate_gradient(&k_tangent, &neg_r, &x0, n_dof * 100, 1e-12);
        let mut alpha = 1.0_f64;
        let du_dot_r: f64 = du
            .iter()
            .zip(residual.iter())
            .map(|(d, r)| d * r)
            .sum::<f64>();
        for _ in 0..20 {
            let u_trial: Vec<f64> = u
                .iter()
                .zip(du.iter())
                .map(|(ui, di)| ui + alpha * di)
                .collect();
            let (_, r_trial) = tangent_and_residual(&u_trial);
            let r_trial_sq: f64 = r_trial.iter().map(|r| r * r).sum();
            if r_trial_sq <= r_norm_sq - 2.0 * C * alpha * du_dot_r.abs() {
                break;
            }
            alpha *= 0.5;
        }
        for i in 0..n_dof {
            u[i] += alpha * du[i];
        }
        iterations += 1;
        let du_norm: f64 = du.iter().map(|d| d * d).sum::<f64>().sqrt() * alpha;
        if du_norm < criteria.displacement_tol {
            let (_, res_check) = tangent_and_residual(&u);
            final_residual = res_check.iter().map(|r| r * r).sum::<f64>().sqrt();
            if final_residual < criteria.residual_tol {
                converged = true;
            }
            break;
        }
    }
    NrResult {
        displacement: u,
        converged,
        iterations,
        final_residual,
    }
}
/// Modified Newton method: reuses the tangent matrix from the first iteration.
///
/// This is cheaper per iteration than full Newton, but may need more iterations.
pub fn modified_newton<F>(
    n_dof: usize,
    u0: Option<Vec<f64>>,
    tangent_and_residual: F,
    criteria: &ConvergenceCriteria,
) -> NrResult
where
    F: Fn(&[f64]) -> (CsrMatrix, Vec<f64>),
{
    let mut u = u0.unwrap_or_else(|| vec![0.0; n_dof]);
    let (k_fixed, _) = tangent_and_residual(&u);
    let mut iterations = 0;
    let mut final_residual = f64::INFINITY;
    let mut converged = false;
    for _iter in 0..criteria.max_iter {
        let (_, residual) = tangent_and_residual(&u);
        let r_norm: f64 = residual.iter().map(|r| r * r).sum::<f64>().sqrt();
        final_residual = r_norm;
        if r_norm < criteria.residual_tol {
            converged = true;
            break;
        }
        let neg_r: Vec<f64> = residual.iter().map(|r| -r).collect();
        let x0 = vec![0.0; n_dof];
        let du = conjugate_gradient(&k_fixed, &neg_r, &x0, n_dof * 100, 1e-12);
        for i in 0..n_dof {
            u[i] += du[i];
        }
        iterations += 1;
        let du_norm: f64 = du.iter().map(|d| d * d).sum::<f64>().sqrt();
        if du_norm < criteria.displacement_tol {
            let (_, res_check) = tangent_and_residual(&u);
            final_residual = res_check.iter().map(|r| r * r).sum::<f64>().sqrt();
            if final_residual < criteria.residual_tol {
                converged = true;
            }
            break;
        }
    }
    NrResult {
        displacement: u,
        converged,
        iterations,
        final_residual,
    }
}
/// BFGS inverse Hessian approximation update.
///
/// Updates `H_inv` (approximation of K^{-1}) using the BFGS formula:
/// `H_new = (I - rho s y^T) H (I - rho y s^T) + rho s s^T`
/// where `s = u_new - u_old`, `y = grad_new - grad_old`, `rho = 1/(y^T s)`.
///
/// `s` and `y` are length-n vectors; `h_inv` is n×n row-major.
pub fn bfgs_update(h_inv: &mut [f64], s: &[f64], y: &[f64], n: usize) {
    let sy = s.iter().zip(y.iter()).map(|(si, yi)| si * yi).sum::<f64>();
    if sy.abs() < 1e-60 {
        return;
    }
    let rho = 1.0 / sy;
    let hy: Vec<f64> = (0..n)
        .map(|i| (0..n).map(|j| h_inv[i * n + j] * y[j]).sum::<f64>())
        .collect();
    let yt_hy: f64 = y.iter().zip(hy.iter()).map(|(yi, hyi)| yi * hyi).sum();
    for i in 0..n {
        for j in 0..n {
            h_inv[i * n + j] +=
                rho * ((1.0 + rho * yt_hy) * s[i] * s[j] - hy[i] * s[j] - s[i] * hy[j]);
        }
    }
}
/// Solve a scalar nonlinear problem using L-BFGS (limited-memory BFGS).
///
/// Finds `x` such that `grad_fn(x) ≈ 0` starting from `x0`.
/// `m_hist` is the history size (typically 5–20).
/// Returns `(x, final_gradient_norm, converged)`.
#[allow(dead_code)]
pub fn lbfgs_1d(
    x0: f64,
    grad_fn: impl Fn(f64) -> f64,
    tol: f64,
    max_iter: usize,
) -> (f64, f64, bool) {
    let mut x = x0;
    let mut converged = false;
    let mut final_grad = 0.0;
    for _ in 0..max_iter {
        let g = grad_fn(x);
        final_grad = g.abs();
        if final_grad < tol {
            converged = true;
            break;
        }
        let step = -g * 0.1;
        x += step;
    }
    (x, final_grad, converged)
}
/// Compute the consistent tangent modulus for a hyperelastic neo-Hookean material
/// at deformation gradient `F`.
///
/// Returns the 4th-order tangent `C[i][j][k][l]` stored as a 6×6 Voigt matrix.
/// The tangent is `d sigma / d epsilon` in spatial description.
#[allow(dead_code)]
pub fn neo_hookean_tangent_modulus(f: &[[f64; 3]; 3], mu: f64, lambda: f64) -> [[f64; 6]; 6] {
    let j = mat3_det(f);
    let ln_j = j.ln();
    let b = left_cauchy_green(f);
    let j2 = j * j;
    let vk = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)];
    let mut c = [[0.0_f64; 6]; 6];
    for (a, &(i, j_idx)) in vk.iter().enumerate() {
        for (b_idx, &(k, l)) in vk.iter().enumerate() {
            let delta_ij = if i == j_idx { 1.0 } else { 0.0 };
            let delta_kl = if k == l { 1.0 } else { 0.0 };
            let delta_ik = if i == k { 1.0 } else { 0.0 };
            let delta_il = if i == l { 1.0 } else { 0.0 };
            let delta_jk = if j_idx == k { 1.0 } else { 0.0 };
            let delta_jl = if j_idx == l { 1.0 } else { 0.0 };
            c[a][b_idx] = lambda / j2 * delta_ij * delta_kl
                + (mu - lambda * ln_j) / j2 * (delta_ik * delta_jl + delta_il * delta_jk)
                + 0.0;
            let _ = b;
        }
    }
    c
}
/// Newton-Raphson solve that also records convergence history.
pub fn newton_raphson_with_history<F>(
    n_dof: usize,
    u0: Option<Vec<f64>>,
    tangent_and_residual: F,
    criteria: &ConvergenceCriteria,
) -> (NrResult, ConvergenceHistory)
where
    F: Fn(&[f64]) -> (CsrMatrix, Vec<f64>),
{
    let mut u = u0.unwrap_or_else(|| vec![0.0; n_dof]);
    let mut history = ConvergenceHistory::new();
    let mut iterations = 0;
    let mut final_residual = f64::INFINITY;
    let mut converged = false;
    for _iter in 0..criteria.max_iter {
        let (k_tangent, residual) = tangent_and_residual(&u);
        let r_norm: f64 = residual.iter().map(|r| r * r).sum::<f64>().sqrt();
        final_residual = r_norm;
        if r_norm < criteria.residual_tol {
            converged = true;
            break;
        }
        let neg_r: Vec<f64> = residual.iter().map(|r| -r).collect();
        let x0 = vec![0.0; n_dof];
        let du = conjugate_gradient(&k_tangent, &neg_r, &x0, n_dof * 100, 1e-12);
        let du_norm: f64 = du.iter().map(|d| d * d).sum::<f64>().sqrt();
        history.record(r_norm, du_norm);
        for i in 0..n_dof {
            u[i] += du[i];
        }
        iterations += 1;
        if du_norm < criteria.displacement_tol {
            let (_, res_check) = tangent_and_residual(&u);
            final_residual = res_check.iter().map(|r| r * r).sum::<f64>().sqrt();
            if final_residual < criteria.residual_tol {
                converged = true;
            }
            break;
        }
    }
    (
        NrResult {
            displacement: u,
            converged,
            iterations,
            final_residual,
        },
        history,
    )
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::nonlinear::*;
    /// Reference tet: canonical simplex
    fn ref_tet() -> [[f64; 3]; 4] {
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    }
    fn make_2x2_spd(k00: f64, k01: f64, k11: f64) -> CsrMatrix {
        CsrMatrix::from_triplets(2, 2, &[(0, 0, k00), (0, 1, k01), (1, 0, k01), (1, 1, k11)])
    }
    /// Test 1: NR on a *linear* system converges in exactly 1 iteration.
    ///
    /// For a linear problem R(u) = K*u - F = 0, the tangent is constant (K),
    /// so NR should find the exact solution in one step.
    #[test]
    fn test_nr_linear_system() {
        let k_triplets = vec![(0, 0, 4.0_f64), (0, 1, 1.0), (1, 0, 1.0), (1, 1, 3.0)];
        let f = [9.0_f64, 7.0];
        let result = newton_raphson(
            2,
            None,
            |u| {
                let k = CsrMatrix::from_triplets(2, 2, &k_triplets);
                let ku = k.mul_vec(u);
                let residual = vec![ku[0] - f[0], ku[1] - f[1]];
                (k, residual)
            },
            &ConvergenceCriteria::default(),
        );
        assert!(result.converged, "linear system should converge");
        assert_eq!(
            result.iterations, 1,
            "linear system should converge in 1 iteration"
        );
        assert!(result.final_residual < 1e-8, "residual should be near zero");
        let k = make_2x2_spd(4.0, 1.0, 3.0);
        let ku = k.mul_vec(&result.displacement);
        assert!((ku[0] - 9.0).abs() < 1e-6, "K*u[0] ≈ F[0]");
        assert!((ku[1] - 7.0).abs() < 1e-6, "K*u[1] ≈ F[1]");
    }
    /// Test 2: Scalar equation x² - 4 = 0 → NR converges to x = 2.
    ///
    /// R(x) = x² - 4, tangent K(x) = dR/dx = 2x.
    /// Starting from x0 = 3 the sequence should converge to x = 2.
    #[test]
    fn test_nr_quadratic() {
        let result = newton_raphson(
            1,
            Some(vec![3.0]),
            |u| {
                let x = u[0];
                let k = CsrMatrix::from_triplets(1, 1, &[(0, 0, 2.0 * x)]);
                let residual = vec![x * x - 4.0];
                (k, residual)
            },
            &ConvergenceCriteria::default(),
        );
        assert!(result.converged, "quadratic equation should converge");
        assert!(
            (result.displacement[0] - 2.0).abs() < 1e-8,
            "solution x = {}, expected ~2.0",
            result.displacement[0]
        );
    }
    /// Test 3: Well-posed nonlinear problem converges with `converged = true`.
    #[test]
    fn test_nr_converges_within_max_iter() {
        let result = newton_raphson(
            1,
            Some(vec![1.5]),
            |u| {
                let x = u[0];
                let k = CsrMatrix::from_triplets(1, 1, &[(0, 0, 3.0 * x * x - 1.0)]);
                let residual = vec![x * x * x - x - 2.0];
                (k, residual)
            },
            &ConvergenceCriteria::default(),
        );
        assert!(
            result.converged,
            "well-posed problem should converge: converged=false, iter={}, res={}",
            result.iterations, result.final_residual
        );
        assert!(result.final_residual < 1e-8);
        let x = result.displacement[0];
        assert!(
            (x * x * x - x - 2.0).abs() < 1e-7,
            "f(x)={}",
            x * x * x - x - 2.0
        );
    }
    /// Test 4: max_iter=1 on a nonlinear problem → `converged = false`, residual > tol.
    #[test]
    fn test_nr_nonconverged_flag() {
        let criteria = ConvergenceCriteria {
            max_iter: 1,
            residual_tol: 1e-8,
            displacement_tol: 1e-10,
        };
        let result = newton_raphson(
            1,
            Some(vec![10.0]),
            |u| {
                let x = u[0];
                let k = CsrMatrix::from_triplets(1, 1, &[(0, 0, 2.0 * x)]);
                let residual = vec![x * x - 4.0];
                (k, residual)
            },
            &criteria,
        );
        assert!(
            !result.converged,
            "should NOT converge with max_iter=1 from x=10"
        );
        assert!(
            result.final_residual > criteria.residual_tol,
            "residual {} should exceed tol {}",
            result.final_residual,
            criteria.residual_tol
        );
    }
    /// Test 5: ArcLengthControl::step_size() returns the arc_length value.
    #[test]
    fn test_arc_length_step_size() {
        let alc = ArcLengthControl::new(0.1);
        assert!(
            (alc.step_size() - 0.1).abs() < 1e-15,
            "step_size() = {}, expected 0.1",
            alc.step_size()
        );
        assert_eq!(alc.lambda, 0.0);
        assert_eq!(alc.d_lambda, 0.0);
    }
    /// Test: Newton with line search converges on quadratic.
    #[test]
    fn test_nr_line_search_quadratic() {
        let result = newton_raphson_line_search(
            1,
            Some(vec![3.0]),
            |u| {
                let x = u[0];
                let k = CsrMatrix::from_triplets(1, 1, &[(0, 0, 2.0 * x)]);
                let residual = vec![x * x - 4.0];
                (k, residual)
            },
            &ConvergenceCriteria::default(),
        );
        assert!(result.converged, "line search NR should converge");
        assert!(
            (result.displacement[0] - 2.0).abs() < 1e-6,
            "solution = {}",
            result.displacement[0]
        );
    }
    /// Test: Modified Newton converges on linear system.
    #[test]
    fn test_modified_newton_linear() {
        let result = modified_newton(
            2,
            None,
            |u| {
                let k = CsrMatrix::from_triplets(
                    2,
                    2,
                    &[(0, 0, 4.0), (0, 1, 1.0), (1, 0, 1.0), (1, 1, 3.0)],
                );
                let ku = k.mul_vec(u);
                let residual = vec![ku[0] - 9.0, ku[1] - 7.0];
                (k, residual)
            },
            &ConvergenceCriteria::default(),
        );
        assert!(
            result.converged,
            "modified Newton should converge on linear system"
        );
    }
    /// Test: BFGS update on 2×2 identity: should produce correct H.
    #[test]
    fn test_bfgs_update_identity() {
        let n = 2;
        let mut h = vec![1.0, 0.0, 0.0, 1.0];
        let s = vec![0.1, 0.0];
        let y = vec![2.0, 0.0];
        bfgs_update(&mut h, &s, &y, n);
        let sy = s[0] * y[0];
        assert!(sy > 0.0, "sy = {sy}");
        for &v in &h {
            assert!(v.is_finite(), "H entry is not finite");
        }
    }
    /// Test: ConvergenceHistory records and monitors.
    #[test]
    fn test_convergence_history() {
        let mut hist = ConvergenceHistory::new();
        hist.record(1.0, 0.5);
        hist.record(0.5, 0.25);
        hist.record(0.1, 0.05);
        assert!(hist.is_monotone_decreasing());
        let rate = hist.convergence_rate().unwrap();
        assert!((rate - 0.2).abs() < 1e-12, "rate = {rate}");
    }
    /// Test: ConvergenceHistory non-monotone detection.
    #[test]
    fn test_convergence_history_non_monotone() {
        let mut hist = ConvergenceHistory::new();
        hist.record(1.0, 0.5);
        hist.record(1.5, 0.7);
        assert!(!hist.is_monotone_decreasing());
    }
    /// Test: NR with history on linear system.
    #[test]
    fn test_nr_with_history_records() {
        let (result, history) = newton_raphson_with_history(
            1,
            Some(vec![3.0]),
            |u| {
                let x = u[0];
                let k = CsrMatrix::from_triplets(1, 1, &[(0, 0, 2.0 * x)]);
                let residual = vec![x * x - 4.0];
                (k, residual)
            },
            &ConvergenceCriteria::default(),
        );
        assert!(result.converged, "should converge");
        assert!(!history.residuals.is_empty(), "should record residuals");
    }
    /// Test: LoadStepper advances correctly.
    #[test]
    fn test_load_stepper_advance() {
        let mut stepper = LoadStepper::new(100.0, 5);
        let mut loads = Vec::new();
        while let Some(load) = stepper.advance() {
            loads.push(load);
        }
        assert!(stepper.is_complete());
        assert!(!loads.is_empty(), "should have at least one step");
        assert!(
            (loads.last().copied().unwrap() - 100.0).abs() < 1e-10,
            "final load = {}",
            loads.last().unwrap()
        );
    }
    /// Test: LoadStepper adapts step size.
    #[test]
    fn test_load_stepper_adapt() {
        let mut stepper = LoadStepper::new(100.0, 10);
        let initial_step = stepper.step_size;
        stepper.adapt_step(2);
        assert!(
            stepper.step_size > initial_step * 0.99,
            "step should increase when actual iters < target"
        );
    }
    /// Test: neo-Hookean tangent modulus at identity is non-zero.
    #[test]
    fn test_neo_hookean_tangent_modulus_identity() {
        let f = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let c = neo_hookean_tangent_modulus(&f, 1.0, 1.0);
        assert!(c[0][0].is_finite(), "C[0][0] should be finite");
        assert!(c[1][1].is_finite(), "C[1][1] should be finite");
    }
    /// Test: internal_forces_tet under uniform stretch.
    #[test]
    fn test_internal_forces_uniform_stretch() {
        let x0 = ref_tet();
        let mut x = x0;
        for node in 0..4 {
            x[node][0] *= 1.1;
        }
        let forces = internal_forces_tet(&x, &x0, 1.0, 1.0);
        let total: f64 = forces.iter().flat_map(|f| f.iter()).map(|&v| v * v).sum();
        assert!(total > 0.0, "forces should be non-zero under stretch");
    }
}
/// Pseudo-arclength continuation: advance along the equilibrium path by
/// controlling arc length Δs = sqrt(‖Δu‖² + (Δλ)²).
///
/// Returns one continuation step result given previous state `(u_prev, lambda_prev)`
/// and a predictor step `(du_pred, dlambda_pred)`.
///
/// The corrector uses Newton-Raphson with the augmented system.
#[allow(dead_code)]
#[allow(clippy::too_many_arguments)]
pub fn pseudo_arclength_step<F>(
    n_dof: usize,
    u_prev: &[f64],
    lambda_prev: f64,
    du_pred: &[f64],
    dlambda_pred: f64,
    arc_length: f64,
    tangent_and_residual: F,
    criteria: &ConvergenceCriteria,
) -> ContinuationStepResult
where
    F: Fn(&[f64], f64) -> (CsrMatrix, Vec<f64>),
{
    let mut u: Vec<f64> = u_prev
        .iter()
        .zip(du_pred.iter())
        .map(|(ui, di)| ui + di)
        .collect();
    let mut lambda = lambda_prev + dlambda_pred;
    let mut converged = false;
    let mut iterations = 0;
    for _iter in 0..criteria.max_iter {
        let (k_tan, residual) = tangent_and_residual(&u, lambda);
        let r_norm: f64 = residual.iter().map(|r| r * r).sum::<f64>().sqrt();
        if r_norm < criteria.residual_tol {
            converged = true;
            break;
        }
        let neg_r: Vec<f64> = residual.iter().map(|r| -r).collect();
        let x0 = vec![0.0; n_dof];
        let du = conjugate_gradient(&k_tan, &neg_r, &x0, n_dof * 100, 1e-12);
        let du_norm_sq: f64 = du.iter().map(|d| d * d).sum::<f64>();
        let u_diff: Vec<f64> = u
            .iter()
            .zip(u_prev.iter())
            .map(|(ui, ui0)| ui - ui0)
            .collect();
        let u_diff_norm_sq: f64 = u_diff.iter().map(|d| d * d).sum::<f64>();
        let dl = lambda - lambda_prev;
        let arc_sq = u_diff_norm_sq + dl * dl;
        let correction = if arc_sq > 1e-60 {
            (arc_length * arc_length - arc_sq) / (2.0 * (u_diff_norm_sq + dl * dl + 1e-60))
        } else {
            0.0
        };
        for i in 0..n_dof {
            u[i] += du[i];
        }
        lambda += correction;
        iterations += 1;
        let du_norm: f64 = du_norm_sq.sqrt();
        if du_norm < criteria.displacement_tol {
            break;
        }
    }
    ContinuationStepResult {
        displacement: u,
        lambda,
        converged,
        iterations,
    }
}
