//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{NeoHookean, OgdenTerm};

/// Multiply two 3x3 matrices: result\[i\]\[j\] = Σ_k a\[i\]\[k\] * b\[k\]\[j\].
pub fn mat3_mul(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut c = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}
/// Transpose a 3x3 matrix.
pub fn mat3_transpose(a: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut t = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            t[i][j] = a[j][i];
        }
    }
    t
}
/// Determinant of a 3x3 matrix (Sarrus' rule).
pub fn mat3_det(a: [[f64; 3]; 3]) -> f64 {
    a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
}
/// Inverse of a 3x3 matrix. Returns `None` if singular (|det| < 1e-14).
pub fn mat3_inverse(a: [[f64; 3]; 3]) -> Option<[[f64; 3]; 3]> {
    let det = mat3_det(a);
    if det.abs() < 1e-14 {
        return None;
    }
    let inv_det = 1.0 / det;
    let mut inv = [[0.0_f64; 3]; 3];
    inv[0][0] = (a[1][1] * a[2][2] - a[1][2] * a[2][1]) * inv_det;
    inv[0][1] = (a[0][2] * a[2][1] - a[0][1] * a[2][2]) * inv_det;
    inv[0][2] = (a[0][1] * a[1][2] - a[0][2] * a[1][1]) * inv_det;
    inv[1][0] = (a[1][2] * a[2][0] - a[1][0] * a[2][2]) * inv_det;
    inv[1][1] = (a[0][0] * a[2][2] - a[0][2] * a[2][0]) * inv_det;
    inv[1][2] = (a[0][2] * a[1][0] - a[0][0] * a[1][2]) * inv_det;
    inv[2][0] = (a[1][0] * a[2][1] - a[1][1] * a[2][0]) * inv_det;
    inv[2][1] = (a[0][1] * a[2][0] - a[0][0] * a[2][1]) * inv_det;
    inv[2][2] = (a[0][0] * a[1][1] - a[0][1] * a[1][0]) * inv_det;
    Some(inv)
}
/// 3x3 identity matrix.
pub fn mat3_identity() -> [[f64; 3]; 3] {
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
}
/// Right Cauchy-Green deformation tensor C = F^T * F.
pub fn right_cauchy_green(f: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    mat3_mul(mat3_transpose(f), f)
}
/// Left Cauchy-Green deformation tensor b = F * F^T.
pub fn left_cauchy_green(f: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    mat3_mul(f, mat3_transpose(f))
}
/// Trace of a 3x3 matrix.
pub(super) fn mat3_trace(a: [[f64; 3]; 3]) -> f64 {
    a[0][0] + a[1][1] + a[2][2]
}
/// Principal invariants of C:
/// - I1 = tr(C)
/// - I2 = (I1² - tr(C²)) / 2
/// - I3 = det(C)
pub fn invariants(c: [[f64; 3]; 3]) -> (f64, f64, f64) {
    let i1 = mat3_trace(c);
    let c2 = mat3_mul(c, c);
    let i2 = (i1 * i1 - mat3_trace(c2)) / 2.0;
    let i3 = mat3_det(c);
    (i1, i2, i3)
}
/// Scale all entries of a 3x3 matrix by a scalar.
pub(super) fn mat3_scale(a: [[f64; 3]; 3], s: f64) -> [[f64; 3]; 3] {
    let mut b = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            b[i][j] = a[i][j] * s;
        }
    }
    b
}
/// Uniaxial Cauchy stress for a sum of Ogden terms.
///
/// σ = Σ μ_i * λ1^α_i
pub fn ogden_cauchy_stress_uniaxial(terms: &[OgdenTerm], lambda1: f64) -> f64 {
    terms.iter().map(|t| t.mu * lambda1.powf(t.alpha)).sum()
}
/// Deviatoric (isochoric) deformation gradient F̄ = J^(-1/3) * F.
pub fn deviatoric_f(f: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let j = mat3_det(f);
    let scale = j.powf(-1.0 / 3.0);
    mat3_scale(f, scale)
}
/// Volumetric strain energy: K/4 * (J² - 1 - 2*ln(J)).
pub fn volumetric_strain_energy(k: f64, j: f64) -> f64 {
    k / 4.0 * (j * j - 1.0 - 2.0 * j.ln())
}
/// Compute the deformation gradient F for a linear tetrahedral element.
///
/// Given reference and current coordinates of 4 nodes, F = dx/dX is
/// computed from the shape function gradients of the linear tetrahedron.
pub fn deformation_gradient_tet(
    ref_coords: &[[f64; 3]; 4],
    cur_coords: &[[f64; 3]; 4],
) -> [[f64; 3]; 3] {
    let mut x_mat = [[0.0_f64; 3]; 3];
    let mut u_mat = [[0.0_f64; 3]; 3];
    for col in 0..3 {
        for row in 0..3 {
            x_mat[row][col] = ref_coords[col + 1][row] - ref_coords[0][row];
            u_mat[row][col] = cur_coords[col + 1][row] - cur_coords[0][row];
        }
    }
    let x_inv = mat3_inverse(x_mat).unwrap_or(mat3_identity());
    mat3_mul(u_mat, x_inv)
}
/// Compute the first Piola-Kirchhoff stress P = dW/dF by central finite differences.
pub(super) fn numerical_pk1(
    f: [[f64; 3]; 3],
    energy: impl Fn([[f64; 3]; 3]) -> f64,
) -> [[f64; 3]; 3] {
    let h = 1e-7_f64;
    let mut p = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let mut f_plus = f;
            let mut f_minus = f;
            f_plus[i][j] += h;
            f_minus[i][j] -= h;
            p[i][j] = (energy(f_plus) - energy(f_minus)) / (2.0 * h);
        }
    }
    p
}
/// Compute the second Piola-Kirchhoff stress S from the first Piola-Kirchhoff P.
///
/// S = F^{-1} P (for symmetric materials: S = F^{-1} dW/dF).
pub fn pk1_to_pk2(f: [[f64; 3]; 3], p: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let f_inv = mat3_inverse(f).unwrap_or(mat3_identity());
    mat3_mul(f_inv, p)
}
/// Compute the Cauchy stress sigma from the first Piola-Kirchhoff P.
///
/// sigma = P F^T / J
pub fn pk1_to_cauchy(f: [[f64; 3]; 3], p: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let j = mat3_det(f);
    let ft = mat3_transpose(f);
    let pft = mat3_mul(p, ft);
    mat3_scale(pft, 1.0 / j)
}
/// Compute the material tangent (4th-order elasticity tensor) numerically.
///
/// Returns a flattened 9x9 matrix `C[iI][jJ]` where `iI = 3*i + I`, `jJ = 3*j + J`.
/// C_{iIjJ} = d²W / (dF_{iI} dF_{jJ})
pub fn material_tangent(f: [[f64; 3]; 3], energy: impl Fn([[f64; 3]; 3]) -> f64) -> Vec<f64> {
    let h = 1e-6_f64;
    let mut tangent = vec![0.0_f64; 81];
    let w0 = energy(f);
    for i in 0..3 {
        for big_i in 0..3 {
            let ii = 3 * i + big_i;
            for j in 0..3 {
                for big_j in 0..3 {
                    let jj = 3 * j + big_j;
                    let mut f_pp = f;
                    let mut f_pm = f;
                    let mut f_mp = f;
                    let mut f_mm = f;
                    f_pp[i][big_i] += h;
                    f_pp[j][big_j] += h;
                    f_pm[i][big_i] += h;
                    f_pm[j][big_j] -= h;
                    f_mp[i][big_i] -= h;
                    f_mp[j][big_j] += h;
                    f_mm[i][big_i] -= h;
                    f_mm[j][big_j] -= h;
                    let _ = w0;
                    tangent[ii * 9 + jj] =
                        (energy(f_pp) - energy(f_pm) - energy(f_mp) + energy(f_mm)) / (4.0 * h * h);
                }
            }
        }
    }
    tangent
}
/// Add a 3x3 matrix `a` to another: `a[i][j] += b[i][j]`.
#[cfg(test)]
pub(super) fn mat3_add(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut c = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] = a[i][j] + b[i][j];
        }
    }
    c
}
/// Frobenius norm of a 3x3 matrix.
#[cfg(test)]
pub(super) fn mat3_frobenius(a: [[f64; 3]; 3]) -> f64 {
    let mut sum = 0.0;
    for row in &a {
        for &val in row {
            sum += val * val;
        }
    }
    sum.sqrt()
}
/// Green-Lagrange strain tensor E = 1/2 * (C - I).
pub fn green_lagrange_strain(f: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let c = right_cauchy_green(f);
    let i = mat3_identity();
    let mut e = [[0.0_f64; 3]; 3];
    for ii in 0..3 {
        for jj in 0..3 {
            e[ii][jj] = 0.5 * (c[ii][jj] - i[ii][jj]);
        }
    }
    e
}
/// Compute the linearised strain (small-strain) tensor from displacement gradient ∇u.
///
/// ε = 1/2 * (∇u + ∇u^T)
pub fn linearised_strain(grad_u: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut eps = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            eps[i][j] = 0.5 * (grad_u[i][j] + grad_u[j][i]);
        }
    }
    eps
}
/// Total Lagrangian internal virtual work integrand: P : δF.
///
/// Returns the scalar P_{iI} δF_{iI}.
pub fn internal_virtual_work(p: [[f64; 3]; 3], df: [[f64; 3]; 3]) -> f64 {
    let mut work = 0.0;
    for i in 0..3 {
        for j in 0..3 {
            work += p[i][j] * df[i][j];
        }
    }
    work
}
/// Push-forward of a 2nd-order material tensor T from reference to current config.
///
/// t_ij = 1/J * F_{iI} T_{IJ} F_{jJ}
pub fn push_forward_tensor(f: [[f64; 3]; 3], t: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let j = mat3_det(f);
    if j.abs() < 1e-14 {
        return t;
    }
    let ft = mat3_transpose(f);
    let fft = mat3_mul(f, mat3_mul(t, ft));
    mat3_scale(fft, 1.0 / j)
}
/// Pull-back of a Cauchy stress tensor σ to the 2nd Piola-Kirchhoff stress S.
///
/// S = J * F^{-1} σ F^{-T}
pub fn pull_back_cauchy_to_pk2(f: [[f64; 3]; 3], sigma: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let j = mat3_det(f);
    let f_inv = mat3_inverse(f).unwrap_or(mat3_identity());
    let f_inv_t = mat3_transpose(f_inv);
    let tmp = mat3_mul(f_inv, sigma);
    mat3_scale(mat3_mul(tmp, f_inv_t), j)
}
/// Deviatoric right Cauchy-Green tensor C̄ = J^{-2/3} C.
pub fn deviatoric_right_cauchy_green(f: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let c = right_cauchy_green(f);
    let j = mat3_det(f);
    let scale = j.powf(-2.0 / 3.0);
    mat3_scale(c, scale)
}
/// Deviatoric part of the strain energy (Neo-Hookean deviatoric).
///
/// W_dev = μ/2 * (Ī1 - 3) where Ī1 = J^{-2/3} * I1
pub fn neo_hookean_deviatoric_energy(mu: f64, f: [[f64; 3]; 3]) -> f64 {
    let c_bar = deviatoric_right_cauchy_green(f);
    let (i1_bar, _, _) = invariants(c_bar);
    mu / 2.0 * (i1_bar - 3.0)
}
/// Deviatoric PK1 stress contribution from Neo-Hookean.
///
/// P_dev = μ * J^{-2/3} * (F - I1/3 * F^{-T})
pub fn neo_hookean_deviatoric_pk1(mu: f64, f: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let j = mat3_det(f);
    let j_neg23 = j.powf(-2.0 / 3.0);
    let c = right_cauchy_green(f);
    let (i1, _, _) = invariants(c);
    let f_inv = mat3_inverse(f).unwrap_or(mat3_identity());
    let f_inv_t = mat3_transpose(f_inv);
    let mut p = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j_idx in 0..3 {
            p[i][j_idx] = mu * j_neg23 * (f[i][j_idx] - i1 / 3.0 * f_inv_t[i][j_idx]);
        }
    }
    p
}
/// B-bar volume-averaged deformation gradient for a single tetrahedral element.
///
/// For the element with integration point deformation gradients F_ip,
/// the B-bar F is F_bar = (J_avg / J_ip)^{1/3} * F_ip.
pub fn fbar_deformation_gradient(f_ip: [[f64; 3]; 3], j_avg: f64) -> [[f64; 3]; 3] {
    let j_ip = mat3_det(f_ip);
    if j_ip.abs() < 1e-14 {
        return f_ip;
    }
    let scale = (j_avg / j_ip).powf(1.0 / 3.0);
    mat3_scale(f_ip, scale)
}
/// F-bar modified strain energy for Neo-Hookean.
///
/// Uses F-bar = (J_avg/J)^{1/3} * F so that the volumetric part uses J_avg.
pub fn neo_hookean_fbar_energy(mu: f64, lambda: f64, f: [[f64; 3]; 3], j_avg: f64) -> f64 {
    let f_bar = fbar_deformation_gradient(f, j_avg);
    let c_bar = right_cauchy_green(f_bar);
    let (i1_bar, _, _) = invariants(c_bar);
    let _j = mat3_det(f);
    let ln_j = j_avg.ln();
    mu / 2.0 * (i1_bar - 3.0) - mu * ln_j + lambda / 2.0 * ln_j * ln_j
}
/// Shape function gradient for a single tetrahedral element at node 0.
///
/// Returns ∇N_0 = -\[N1+N2+N3\] gradient, which equals minus the sum of
/// gradients of nodes 1,2,3 (partition of unity).
pub fn tet_shape_gradient_node0(ref_coords: &[[f64; 3]; 4]) -> [f64; 3] {
    let mut grad = [0.0_f64; 3];
    for node in 1..4 {
        for dim in 0..3 {
            grad[dim] -= ref_coords[node][dim] - ref_coords[0][dim];
        }
    }
    grad
}
/// Check Drucker-Prager stability condition for a hyperelastic material.
///
/// The material is stable at F if dP:dF ≥ 0 for all small perturbations dF.
/// This checks a random set of 6 axis-aligned perturbations.
///
/// Returns `true` if stable (all dP:dF ≥ -tol).
pub fn check_drucker_stability(
    f: [[f64; 3]; 3],
    energy: &impl Fn([[f64; 3]; 3]) -> f64,
    tolerance: f64,
) -> bool {
    let h = 1e-5_f64;
    let p0 = numerical_pk1(f, energy);
    for i in 0..3 {
        for j in 0..3 {
            let mut f_pert = f;
            f_pert[i][j] += h;
            let p_pert = numerical_pk1(f_pert, energy);
            let dp = p_pert[i][j] - p0[i][j];
            let df = h;
            if dp * df < -tolerance {
                return false;
            }
        }
    }
    true
}
/// Compute the Legendre-Hadamard (strong ellipticity) condition.
///
/// The material is strongly elliptic if:
/// (m⊗n) : A : (m⊗n) > 0
/// for all non-zero vectors m (material) and n (spatial).
///
/// Returns the minimum acoustic tensor eigenvalue estimate over a set
/// of direction pairs. A positive value indicates strong ellipticity.
pub fn legendre_hadamard_check(f: [[f64; 3]; 3], tangent: &[f64]) -> f64 {
    let dirs: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let mut min_val = f64::INFINITY;
    let _ = f;
    for n in &dirs {
        for m in &dirs {
            let mut q = [[0.0_f64; 3]; 3];
            for (i, q_row) in q.iter_mut().enumerate() {
                for (k, q_ik) in q_row.iter_mut().enumerate() {
                    for big_i in 0..3 {
                        for big_j in 0..3 {
                            let idx_a = (i * 3 + big_i) * 9 + k * 3 + big_j;
                            if idx_a < tangent.len() {
                                *q_ik += tangent[idx_a] * n[big_i] * n[big_j];
                            }
                        }
                    }
                }
            }
            let mut val = 0.0;
            for i in 0..3 {
                for k in 0..3 {
                    val += m[i] * q[i][k] * m[k];
                }
            }
            min_val = min_val.min(val);
        }
    }
    min_val
}
/// Energy-based stability indicator for a hyperelastic material.
///
/// Returns the ratio W(F_perturbed) / W(F) for a small isochoric perturbation.
/// Values > 1 indicate the material resists isochoric distortions (stable).
pub fn isochoric_stability_ratio(
    f: [[f64; 3]; 3],
    energy: &impl Fn([[f64; 3]; 3]) -> f64,
    perturbation: f64,
) -> f64 {
    let w0 = energy(f);
    if w0.abs() < 1e-30 {
        return 1.0;
    }
    let mut f_pert = f;
    f_pert[0][1] += perturbation;
    let w_pert = energy(f_pert);
    w_pert / w0
}
/// Compute the consistent (algorithmic) tangent moduli as a 9×9 matrix.
///
/// This is equivalent to `material_tangent` but specifically designed for
/// use with hyperelastic models that implement the `NeoHookean`-like interface.
///
/// Returns `C[iI][jJ]` where the flattened index is `iI = 3*i + I`.
pub fn consistent_tangent_4th_order(nh: &NeoHookean, f: [[f64; 3]; 3]) -> Vec<f64> {
    material_tangent(f, |ff| nh.strain_energy(ff))
}
/// Compute the consistent Lagrangian elasticity tensor as a 6×6 Voigt matrix.
///
/// Maps 2nd Piola-Kirchhoff stress increments to Green-Lagrange strain increments.
/// The Voigt ordering is: \[11, 22, 33, 12, 23, 13\].
pub fn consistent_material_tangent_voigt(nh: &NeoHookean, f: [[f64; 3]; 3]) -> [[f64; 6]; 6] {
    let voigt_map = [(0, 0), (1, 1), (2, 2), (0, 1), (1, 2), (0, 2)];
    let h = 1e-6_f64;
    let c = right_cauchy_green(f);
    let mut cmat = [[0.0_f64; 6]; 6];
    for (beta, &(b_i, b_j)) in voigt_map.iter().enumerate() {
        let mut c_plus = c;
        let mut c_minus = c;
        c_plus[b_i][b_j] += h;
        c_plus[b_j][b_i] += h;
        c_minus[b_i][b_j] -= h;
        c_minus[b_j][b_i] -= h;
        let s_plus = nh_pk2_from_c(nh, c_plus);
        let s_minus = nh_pk2_from_c(nh, c_minus);
        for (alpha, &(a_i, a_j)) in voigt_map.iter().enumerate() {
            let ds = (s_plus[a_i][a_j] - s_minus[a_i][a_j]) / (2.0 * h);
            cmat[alpha][beta] = ds;
        }
    }
    cmat
}
/// Helper: compute PK2 stress for Neo-Hookean given C (via F = sqrt(C) approximation).
pub(super) fn nh_pk2_from_c(nh: &NeoHookean, c: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let j2 = mat3_det(c).max(1e-14);
    let j = j2.sqrt();
    let ln_j = j.ln();
    let c_inv = mat3_inverse(c).unwrap_or(mat3_identity());
    let mut s = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for k in 0..3 {
            let delta = if i == k { 1.0 } else { 0.0 };
            s[i][k] = nh.mu * (delta - c_inv[i][k]) + nh.lambda * ln_j * c_inv[i][k];
        }
    }
    s
}
/// Numerical second Piola-Kirchhoff stress by differentiating W w.r.t. C.
///
/// `S = 2 * dW/dC` using central differences.
pub fn numerical_pk2_from_energy(
    c: [[f64; 3]; 3],
    energy_fn: impl Fn([[f64; 3]; 3]) -> f64,
) -> [[f64; 3]; 3] {
    let h = 1e-7_f64;
    let mut s = [[0.0_f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let mut c_plus = c;
            let mut c_minus = c;
            c_plus[i][j] += h;
            c_plus[j][i] += h;
            c_minus[i][j] -= h;
            c_minus[j][i] -= h;
            s[i][j] = (energy_fn(c_plus) - energy_fn(c_minus)) / (2.0 * h);
        }
    }
    s
}
