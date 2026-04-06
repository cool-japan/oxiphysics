//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use super::operations::*;
use super::types::*;

/// Right polar decomposition of a non-singular 3×3 tensor F = R U.
///
/// Uses iterative algorithm: R_new = (R + R^{-T})/2 normalised.
/// Returns `(R, U)` where R is proper orthogonal and U is symmetric positive-definite.
#[allow(dead_code)]
pub fn polar_decompose_right(f: &Tensor2) -> Option<(Tensor2, Tensor2)> {
    let mut r = Tensor2 { data: f.data };
    for _ in 0..200 {
        let r_inv_t = r.inverse()?.transpose();
        let r_new_data: [[f64; 3]; 3] = {
            let mut d = [[0.0f64; 3]; 3];
            for i in 0..3 {
                for j in 0..3 {
                    d[i][j] = 0.5 * (r.data[i][j] + r_inv_t.data[i][j]);
                }
            }
            d
        };
        let r_new = Tensor2 { data: r_new_data };
        let diff_norm: f64 = (0..3)
            .flat_map(|i| (0..3).map(move |j| (r_new.data[i][j] - r.data[i][j]).powi(2)))
            .sum::<f64>()
            .sqrt();
        r = r_new;
        if diff_norm < 1e-12 {
            break;
        }
    }
    let rt = r.transpose();
    let u = rt.dot(f);
    Some((r, u))
}
/// Approximate matrix logarithm of a tensor close to identity: ln(I + X) ≈ X - X²/2 + X³/3.
///
/// Only accurate for ||X|| < 0.5 (small deformation regime).
#[allow(dead_code)]
pub fn log_tensor_approx(f: &Tensor2) -> Tensor2 {
    let id = Tensor2::identity();
    let x = f.sub(&id);
    let x2 = x.dot(&x);
    let x3 = x2.dot(&x);
    let term1 = Tensor2 { data: x.data };
    let term2 = x2.scale(0.5);
    let term3 = x3.scale(1.0 / 3.0);
    term1.sub(&term2).add(&term3)
}
/// Approximate matrix exponential via Padé (2,2) approximation.
///
/// exp(A) ≈ (I - A/2 + A²/12)^{-1} (I + A/2 + A²/12).
#[allow(dead_code)]
pub fn exp_tensor_pade(a: &Tensor2) -> Option<Tensor2> {
    let id = Tensor2::identity();
    let a2 = a.dot(a);
    let a_half = a.scale(0.5);
    let a2_12 = a2.scale(1.0 / 12.0);
    let p = id.add(&a_half).add(&a2_12);
    let q = id.sub(&a_half).add(&a2_12);
    let q_inv = q.inverse()?;
    Some(q_inv.dot(&p))
}
/// Einstein summation over a pair of DenseTensors with explicit free and contracted indices.
///
/// `free_a`: list of modes of `a` that appear in the output.
/// `free_b`: list of modes of `b` that appear in the output.
/// `contract_a`: list of modes of `a` to contract.
/// `contract_b`: list of modes of `b` to contract (must pair with `contract_a`).
///
/// Output shape = \[a.shape\[free_a\[0\]\], a.shape\[free_a\[1\]\], ..., b.shape\[free_b\[0\]\], ...\].
#[allow(dead_code)]
#[allow(clippy::too_many_arguments)]
pub fn general_einsum(
    a: &DenseTensor,
    free_a: &[usize],
    contract_a: &[usize],
    b: &DenseTensor,
    free_b: &[usize],
    contract_b: &[usize],
) -> DenseTensor {
    assert_eq!(contract_a.len(), contract_b.len());
    for (&ca, &cb) in contract_a.iter().zip(contract_b.iter()) {
        assert_eq!(a.shape[ca], b.shape[cb], "contracted dimensions must match");
    }
    let mut out_shape: Vec<usize> = free_a.iter().map(|&m| a.shape[m]).collect();
    out_shape.extend(free_b.iter().map(|&m| b.shape[m]));
    let n_out: usize = if out_shape.is_empty() {
        1
    } else {
        out_shape.iter().product()
    };
    let mut out_data = vec![0.0f64; n_out];
    let out_rank = out_shape.len();
    let mut out_strides = vec![1usize; out_rank.max(1)];
    for k in (0..out_rank.saturating_sub(1)).rev() {
        out_strides[k] = out_strides[k + 1] * out_shape[k + 1];
    }
    let strides_a = a.strides();
    let strides_b = b.strides();
    for fa in 0..a.data.len() {
        let mut tmp = fa;
        let mut midx_a = vec![0usize; a.shape.len()];
        for k in 0..a.shape.len() {
            midx_a[k] = tmp / strides_a[k];
            tmp %= strides_a[k];
        }
        for fb in 0..b.data.len() {
            let mut tmp2 = fb;
            let mut midx_b = vec![0usize; b.shape.len()];
            for k in 0..b.shape.len() {
                midx_b[k] = tmp2 / strides_b[k];
                tmp2 %= strides_b[k];
            }
            let contracted = contract_a
                .iter()
                .zip(contract_b.iter())
                .all(|(&ca, &cb)| midx_a[ca] == midx_b[cb]);
            if !contracted {
                continue;
            }
            let mut oidx = vec![0usize; out_rank];
            let mut oi = 0;
            for &ma in free_a.iter() {
                oidx[oi] = midx_a[ma];
                oi += 1;
            }
            for &mb in free_b.iter() {
                oidx[oi] = midx_b[mb];
                oi += 1;
            }
            let fo: usize = if out_rank == 0 {
                0
            } else {
                oidx.iter()
                    .zip(out_strides.iter())
                    .map(|(&i, &s)| i * s)
                    .sum()
            };
            out_data[fo] += a.data[fa] * b.data[fb];
        }
    }
    DenseTensor {
        shape: out_shape,
        data: out_data,
    }
}
/// Reconstruct a DenseTensor from a CP decomposition.
#[allow(dead_code)]
pub fn cp_reconstruct(cp: &CpDecomposition) -> DenseTensor {
    let r = cp.lambdas.len();
    let n0 = cp.a.len();
    let n1 = cp.b.len();
    let n2 = cp.c.len();
    let mut out = DenseTensor::zeros(&[n0, n1, n2]);
    for rr in 0..r {
        let lam = cp.lambdas[rr];
        for i in 0..n0 {
            for j in 0..n1 {
                for k in 0..n2 {
                    let v = out.get(&[i, j, k]) + lam * cp.a[i][rr] * cp.b[j][rr] * cp.c[k][rr];
                    out.set(&[i, j, k], v);
                }
            }
        }
    }
    out
}
/// Relative reconstruction error: ||T - T_approx||_F / ||T||_F.
#[allow(dead_code)]
pub fn cp_relative_error(original: &DenseTensor, cp: &CpDecomposition) -> f64 {
    let recon = cp_reconstruct(cp);
    let diff_norm = original.sub_tensor(&recon).frobenius_norm();
    let orig_norm = original.frobenius_norm();
    if orig_norm < 1e-30 {
        return diff_norm;
    }
    diff_norm / orig_norm
}
/// Compute the 4th-order stiffness tensor for Neo-Hookean material at deformation F.
///
/// C_ijkl = lambda delta_ij delta_kl + mu (delta_ik delta_jl + delta_il delta_jk - delta_ij delta_kl/3)
/// (linearised tangent at identity).
#[allow(dead_code)]
pub fn neo_hookean_stiffness(lambda: f64, mu: f64) -> Tensor4 {
    Tensor4::isotropic(lambda, mu)
}
/// Compute Cauchy stress from linear elasticity: sigma = C : epsilon.
#[allow(dead_code)]
pub fn cauchy_stress_linear(c: &Tensor4, epsilon: &Tensor2) -> Tensor2 {
    c.double_contract_2(epsilon)
}
/// Compute the engineering Young's modulus E and Poisson's ratio nu from Lame parameters.
#[allow(dead_code)]
pub fn lame_to_young_poisson(lambda: f64, mu: f64) -> (f64, f64) {
    let e = mu * (3.0 * lambda + 2.0 * mu) / (lambda + mu);
    let nu = lambda / (2.0 * (lambda + mu));
    (e, nu)
}
/// Compute Lame parameters from Young's modulus E and Poisson's ratio nu.
#[allow(dead_code)]
pub fn young_poisson_to_lame(e: f64, nu: f64) -> (f64, f64) {
    let lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let mu = e / (2.0 * (1.0 + nu));
    (lambda, mu)
}
/// Compute the bulk modulus K from Lame parameters.
#[allow(dead_code)]
pub fn bulk_modulus(lambda: f64, mu: f64) -> f64 {
    lambda + 2.0 * mu / 3.0
}
/// Compute the shear modulus (= mu).
#[allow(dead_code)]
pub fn shear_modulus(mu: f64) -> f64 {
    mu
}
/// Build a transversely isotropic stiffness tensor given five elastic constants.
///
/// Symmetry axis = e_3 (index 2).
/// C_1111 = C_2222 = c11, C_3333 = c33,
/// C_1122 = c12, C_1133 = C_2233 = c13,
/// C_2323 = C_1313 = c44, C_1212 = (c11-c12)/2.
#[allow(dead_code)]
#[allow(clippy::too_many_arguments)]
pub fn transversely_isotropic_stiffness(
    c11: f64,
    c33: f64,
    c12: f64,
    c13: f64,
    c44: f64,
) -> Tensor4 {
    let c66 = 0.5 * (c11 - c12);
    let mut t = Tensor4::zero();
    let mut set = |i: usize, j: usize, k: usize, l: usize, v: f64| {
        t.data[i][j][k][l] = v;
        t.data[j][i][k][l] = v;
        t.data[i][j][l][k] = v;
        t.data[j][i][l][k] = v;
        t.data[k][l][i][j] = v;
        t.data[l][k][i][j] = v;
        t.data[k][l][j][i] = v;
        t.data[l][k][j][i] = v;
    };
    set(0, 0, 0, 0, c11);
    set(1, 1, 1, 1, c11);
    set(2, 2, 2, 2, c33);
    set(0, 0, 1, 1, c12);
    set(0, 0, 2, 2, c13);
    set(1, 1, 2, 2, c13);
    set(1, 2, 1, 2, c44);
    set(0, 2, 0, 2, c44);
    set(0, 1, 0, 1, c66);
    let _ = c66;
    t
}
/// Harmonic mean of two stiffness tensors (Reuss bound):
/// C_Reuss = (C_a^{-1} + C_b^{-1})^{-1} / 2  (in Voigt/Kelvin form).
#[allow(dead_code)]
pub fn reuss_average(c_a: &Tensor4, c_b: &Tensor4) -> Option<Tensor4> {
    let ma = KelvinTensor::from_tensor4(c_a);
    let mb = KelvinTensor::from_tensor4(c_b);
    let sa = invert_6x6(&ma)?;
    let sb = invert_6x6(&mb)?;
    let mut sc = [[0.0f64; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            sc[i][j] = 0.5 * (sa[i][j] + sb[i][j]);
        }
    }
    let mc = invert_6x6(&sc)?;
    Some(KelvinTensor::to_tensor4(&mc))
}
/// Arithmetic mean of two stiffness tensors (Voigt bound):
/// C_Voigt = (C_a + C_b) / 2.
#[allow(dead_code)]
pub fn voigt_average(c_a: &Tensor4, c_b: &Tensor4) -> Tensor4 {
    let a = c_a.add(c_b);
    a.scale(0.5)
}
/// Hill average: arithmetic mean of Voigt and Reuss bounds.
#[allow(dead_code)]
pub fn hill_average(c_a: &Tensor4, c_b: &Tensor4) -> Option<Tensor4> {
    let cv = voigt_average(c_a, c_b);
    let cr = reuss_average(c_a, c_b)?;
    Some(cv.add(&cr).scale(0.5))
}
#[cfg(test)]
mod tests_extended {
    use super::*;
    fn approx(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-10
    }
    #[test]
    fn test_bilinear_form_identity() {
        let id = Tensor2::identity();
        let v = [1.0, 2.0, 3.0];
        let w = [4.0, 5.0, 6.0];
        let bf = bilinear_form(&id, &v, &w);
        assert!(approx(bf, 1.0 * 4.0 + 2.0 * 5.0 + 3.0 * 6.0));
    }
    #[test]
    fn test_quadratic_form_identity() {
        let id = Tensor2::identity();
        let v = [3.0, 4.0, 0.0];
        let q = quadratic_form(&id, &v);
        assert!(approx(q, 25.0));
    }
    #[test]
    fn test_quadratic_form_zero_tensor() {
        let z = Tensor2::zero();
        let v = [1.0, 2.0, 3.0];
        assert!(approx(quadratic_form(&z, &v), 0.0));
    }
    #[test]
    fn test_bilinear_form_diagonal() {
        let mut d = Tensor2::zero();
        d.data[0][0] = 2.0;
        d.data[1][1] = 3.0;
        d.data[2][2] = 5.0;
        let v = [1.0, 1.0, 1.0];
        let w = [1.0, 1.0, 1.0];
        assert!(approx(bilinear_form(&d, &v, &w), 10.0));
    }
    #[test]
    fn test_vec_outer_rank() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        let t = vec_outer(&a, &b);
        assert!(approx(t.data[0][1], 1.0));
        assert!(approx(t.data[1][0], 0.0));
    }
    #[test]
    fn test_vec_outer_self_is_symmetric_for_same_vec() {
        let v = [1.0, 2.0, 3.0];
        let t = vec_outer(&v, &v);
        assert!(t.is_symmetric(1e-12));
    }
    #[test]
    fn test_tensor_outer_shape() {
        let a = DenseTensor::zeros(&[2, 3]);
        let b = DenseTensor::zeros(&[4, 5]);
        let out = tensor_outer(&a, &b);
        assert_eq!(out.shape, vec![2, 3, 4, 5]);
    }
    #[test]
    fn test_tensor_outer_values() {
        let mut a = DenseTensor::zeros(&[2]);
        a.set(&[0], 2.0);
        a.set(&[1], 3.0);
        let mut b = DenseTensor::zeros(&[2]);
        b.set(&[0], 5.0);
        b.set(&[1], 7.0);
        let out = tensor_outer(&a, &b);
        assert!(approx(out.get(&[0, 0]), 10.0));
        assert!(approx(out.get(&[1, 1]), 21.0));
        assert!(approx(out.get(&[0, 1]), 14.0));
    }
    #[test]
    fn test_symmetrize_tensor_rank2() {
        let mut t = DenseTensor::zeros(&[3, 3]);
        t.set(&[0, 1], 2.0);
        t.set(&[1, 0], 4.0);
        let s = symmetrize_tensor(&t);
        assert!(approx(s.get(&[0, 1]), 3.0));
        assert!(approx(s.get(&[1, 0]), 3.0));
    }
    #[test]
    fn test_symmetrize_tensor_already_sym() {
        let mut t = DenseTensor::zeros(&[2, 2]);
        t.set(&[0, 0], 1.0);
        t.set(&[0, 1], 3.0);
        t.set(&[1, 0], 3.0);
        t.set(&[1, 1], 2.0);
        let s = symmetrize_tensor(&t);
        assert!(approx(s.get(&[0, 1]), 3.0));
    }
    #[test]
    fn test_dense_tensor_scale() {
        let t = DenseTensor::from_data(&[2, 2], vec![1.0, 2.0, 3.0, 4.0]);
        let s = t.scale(2.0);
        assert!(approx(s.get(&[0, 0]), 2.0));
        assert!(approx(s.get(&[1, 1]), 8.0));
    }
    #[test]
    fn test_dense_tensor_add_sub() {
        let a = DenseTensor::from_data(&[3], vec![1.0, 2.0, 3.0]);
        let b = DenseTensor::from_data(&[3], vec![4.0, 5.0, 6.0]);
        let c = a.add_tensor(&b);
        assert!(approx(c.get(&[2]), 9.0));
        let d = c.sub_tensor(&a);
        assert!(approx(d.get(&[0]), 4.0));
    }
    #[test]
    fn test_dense_tensor_sum() {
        let t = DenseTensor::from_data(&[4], vec![1.0, 2.0, 3.0, 4.0]);
        assert!(approx(t.sum(), 10.0));
    }
    #[test]
    fn test_dense_tensor_max_abs() {
        let t = DenseTensor::from_data(&[3], vec![-5.0, 3.0, 1.0]);
        assert!(approx(t.max_abs(), 5.0));
    }
    #[test]
    fn test_fourth_order_symmetric_identity_symmetry() {
        let is4 = fourth_order_symmetric_identity();
        assert!(is4.has_minor_symmetry(1e-12));
        assert!(is4.has_major_symmetry(1e-12));
    }
    #[test]
    fn test_fourth_order_skew_identity_antisymmetry() {
        let ia4 = fourth_order_skew_identity();
        let result = ia4.double_contract_4(&ia4);
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        assert!(
                            (result.data[i][j][k][l] - ia4.data[i][j][k][l]).abs() < 1e-10,
                            "I^A :: I^A != I^A at [{i}][{j}][{k}][{l}]"
                        );
                    }
                }
            }
        }
    }
    #[test]
    fn test_fourth_order_is_plus_ia_equals_identity() {
        let is4 = fourth_order_symmetric_identity();
        let ia4 = fourth_order_skew_identity();
        let id4 = Tensor4::identity();
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        assert!(
                            (is4.data[i][j][k][l] + ia4.data[i][j][k][l] - id4.data[i][j][k][l])
                                .abs()
                                < 1e-12
                        );
                    }
                }
            }
        }
    }
    #[test]
    fn test_elasticity_stress_isotropic() {
        let c = Tensor4::isotropic(1.0, 1.0);
        let eps = Tensor2::identity().scale(1.0 / 3.0);
        let sigma = elasticity_stress(&c, &eps);
        let s_diag = sigma.data[0][0];
        let expected = 1.0 * 1.0 + 2.0 * 1.0 / 3.0;
        assert!(
            (s_diag - expected).abs() < 1e-10,
            "s_diag={s_diag}, expected={expected}"
        );
    }
    #[test]
    fn test_lame_to_young_poisson_round_trip() {
        let (lambda, mu) = (1.0e10, 1.5e10);
        let (e, nu) = lame_to_young_poisson(lambda, mu);
        let (lambda2, mu2) = young_poisson_to_lame(e, nu);
        assert!((lambda - lambda2).abs() / lambda < 1e-10);
        assert!((mu - mu2).abs() / mu < 1e-10);
    }
    #[test]
    fn test_young_poisson_steel() {
        let (lambda, mu) = young_poisson_to_lame(200e9, 0.3);
        assert!(lambda > 0.0 && mu > 0.0);
        let k = bulk_modulus(lambda, mu);
        let k_expected = 200e9 / (3.0 * (1.0 - 2.0 * 0.3));
        assert!((k - k_expected).abs() / k_expected < 1e-8);
    }
    #[test]
    fn test_invert_6x6_identity() {
        let id6: [[f64; 6]; 6] = {
            let mut m = [[0.0f64; 6]; 6];
            for i in 0..6 {
                m[i][i] = 1.0;
            }
            m
        };
        let inv = invert_6x6(&id6).unwrap();
        for i in 0..6 {
            for j in 0..6 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((inv[i][j] - expected).abs() < 1e-10);
            }
        }
    }
    #[test]
    fn test_invert_6x6_diagonal() {
        let mut m = [[0.0f64; 6]; 6];
        for i in 0..6 {
            m[i][i] = (i + 1) as f64;
        }
        let inv = invert_6x6(&m).unwrap();
        for i in 0..6 {
            assert!((inv[i][i] - 1.0 / (i + 1) as f64).abs() < 1e-10);
        }
    }
    #[test]
    fn test_compliance_round_trip_isotropic() {
        let c = Tensor4::isotropic(1.0e10, 1.5e10);
        let s_mat = compliance_from_stiffness(&c).unwrap();
        let c2 = KelvinTensor::to_tensor4(&s_mat);
        let c_kelvin = KelvinTensor::from_tensor4(&c);
        for i in 0..6 {
            let mut row_sum = 0.0f64;
            for k in 0..6 {
                row_sum += c_kelvin[i][k] * s_mat[k][i];
            }
            assert!((row_sum - 1.0).abs() < 1e-6, "diagonal[{i}]={row_sum}");
        }
        let _ = c2;
    }
    #[test]
    fn test_polar_decompose_identity() {
        let f = Tensor2::identity();
        let (r, u) = polar_decompose_right(&f).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((r.data[i][j] - expected).abs() < 1e-8);
                assert!((u.data[i][j] - expected).abs() < 1e-8);
            }
        }
    }
    #[test]
    fn test_polar_decompose_r_is_orthogonal() {
        let f = Tensor2::new([[1.2, 0.3, 0.0], [0.1, 0.9, 0.0], [0.0, 0.0, 1.0]]);
        let (r, _u) = polar_decompose_right(&f).unwrap();
        let rrt = r.dot(&r.transpose());
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (rrt.data[i][j] - expected).abs() < 1e-7,
                    "R R^T != I at [{i}][{j}]: got {}",
                    rrt.data[i][j]
                );
            }
        }
    }
    #[test]
    fn test_log_tensor_approx_near_identity() {
        let epsilon = 0.01;
        let mut a = Tensor2::zero();
        a.data[0][1] = 1.0;
        a.data[1][0] = 1.0;
        let fa = Tensor2::identity().add(&a.scale(epsilon));
        let log_fa = log_tensor_approx(&fa);
        assert!(
            (log_fa.data[0][1] - epsilon).abs() < 1e-4,
            "log_fa[0][1]={}",
            log_fa.data[0][1]
        );
    }
    #[test]
    fn test_exp_tensor_pade_identity() {
        let z = Tensor2::zero();
        let e = exp_tensor_pade(&z).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (e.data[i][j] - expected).abs() < 1e-10,
                    "exp(0) != I at [{i}][{j}]"
                );
            }
        }
    }
    #[test]
    fn test_eshelby_sphere_symmetry() {
        let s = eshelby_sphere(0.3);
        assert!(s.has_minor_symmetry(1e-12));
    }
    #[test]
    fn test_eshelby_sphere_nu_zero() {
        let s = eshelby_sphere(0.0);
        assert!((s.data[0][0][0][0] - 7.0 / 15.0).abs() < 1e-10);
        assert!((s.data[0][0][1][1] - (-1.0 / 15.0)).abs() < 1e-10);
        assert!((s.data[0][1][0][1] - 4.0 / 15.0).abs() < 1e-10);
    }
    #[test]
    fn test_transversely_isotropic_major_symmetry() {
        let c = transversely_isotropic_stiffness(200.0, 180.0, 60.0, 65.0, 75.0);
        assert!(c.has_major_symmetry(1e-10));
    }
    #[test]
    fn test_transversely_isotropic_minor_symmetry() {
        let c = transversely_isotropic_stiffness(200.0, 180.0, 60.0, 65.0, 75.0);
        assert!(c.has_minor_symmetry(1e-10));
    }
    #[test]
    fn test_voigt_average_is_arithmetic_mean() {
        let c1 = Tensor4::isotropic(1.0, 1.0);
        let c2 = Tensor4::isotropic(3.0, 3.0);
        let cv = voigt_average(&c1, &c2);
        let expected = Tensor4::isotropic(2.0, 2.0);
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        assert!((cv.data[i][j][k][l] - expected.data[i][j][k][l]).abs() < 1e-10);
                    }
                }
            }
        }
    }
    #[test]
    fn test_reuss_average_symmetric_matrices() {
        let c1 = Tensor4::isotropic(1.0, 1.0);
        let c2 = Tensor4::isotropic(1.0, 1.0);
        let cr = reuss_average(&c1, &c2).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        assert!(
                            (cr.data[i][j][k][l] - c1.data[i][j][k][l]).abs() < 1e-4,
                            "[{i}][{j}][{k}][{l}]: cr={}, c1={}",
                            cr.data[i][j][k][l],
                            c1.data[i][j][k][l]
                        );
                    }
                }
            }
        }
    }
    #[test]
    fn test_cp_reconstruct_zero_lambdas() {
        let cp = CpDecomposition {
            lambdas: vec![0.0],
            a: vec![vec![1.0], vec![0.0], vec![0.0]],
            b: vec![vec![1.0], vec![0.0], vec![0.0]],
            c: vec![vec![1.0], vec![0.0], vec![0.0]],
        };
        let recon = cp_reconstruct(&cp);
        assert!(recon.frobenius_norm() < 1e-12);
    }
    #[test]
    fn test_cp_relative_error_exact_reconstruction() {
        let t = DenseTensor::from_data(&[2, 2, 2], vec![1.0, 2.0, 2.0, 4.0, 2.0, 4.0, 4.0, 8.0]);
        let cp = cp_als(&t, 1, 50, 1e-10);
        let err = cp_relative_error(&t, &cp);
        assert!(err < 0.01, "relative error too large: {err}");
    }
    #[test]
    fn test_tucker_reconstruction_error_full_rank() {
        let shape = [3, 3, 3];
        let n: usize = shape.iter().product();
        let data: Vec<f64> = (0..n).map(|x| x as f64).collect();
        let t = DenseTensor::from_data(&shape, data);
        let td = tucker_hosvd(&t, shape);
        let err = td.relative_error(&t);
        assert!(err < 1e-10, "full-rank Tucker error: {err}");
    }
    #[test]
    fn test_tt_to_dense_shape() {
        let t = DenseTensor::from_data(&[2, 3], (0..6).map(|x| x as f64).collect());
        let tt = TensorTrain::from_dense(&t, 6, 1e-14);
        let recon = tt.to_dense();
        assert_eq!(recon.shape, t.shape);
    }
    #[test]
    fn test_tt_scale() {
        let t = DenseTensor::from_data(&[2, 2], vec![1.0, 2.0, 3.0, 4.0]);
        let tt = TensorTrain::from_dense(&t, 4, 1e-14);
        let tt2 = tt.scale(2.0);
        let recon = tt2.to_dense();
        assert!((recon.get(&[0, 0]) - 2.0).abs() < 1e-6);
    }
    #[test]
    fn test_tt_dot_product_self() {
        let t = DenseTensor::from_data(&[2, 2], vec![1.0, 2.0, 3.0, 4.0]);
        let tt = TensorTrain::from_dense(&t, 4, 1e-14);
        let dot = tt.dot_product(&tt);
        let expected = t.data.iter().map(|&x| x * x).sum::<f64>();
        assert!((dot - expected).abs() / expected < 1e-6);
    }
    #[test]
    fn test_general_einsum_matmul() {
        let a = DenseTensor::from_data(&[2, 2], vec![1.0, 2.0, 3.0, 4.0]);
        let b = DenseTensor::from_data(&[2, 2], vec![5.0, 6.0, 7.0, 8.0]);
        let c = general_einsum(&a, &[0], &[1], &b, &[1], &[0]);
        assert!((c.get(&[0, 0]) - 19.0).abs() < 1e-10);
        assert!((c.get(&[0, 1]) - 22.0).abs() < 1e-10);
        assert!((c.get(&[1, 0]) - 43.0).abs() < 1e-10);
        assert!((c.get(&[1, 1]) - 50.0).abs() < 1e-10);
    }
    #[test]
    fn test_general_einsum_inner_product() {
        let a = DenseTensor::from_data(&[3], vec![1.0, 2.0, 3.0]);
        let b = DenseTensor::from_data(&[3], vec![4.0, 5.0, 6.0]);
        let c = general_einsum(&a, &[], &[0], &b, &[], &[0]);
        assert!((c.get(&[0]) - 32.0).abs() < 1e-10);
    }
    #[test]
    fn test_general_einsum_outer_product() {
        let a = DenseTensor::from_data(&[2], vec![2.0, 3.0]);
        let b = DenseTensor::from_data(&[2], vec![5.0, 7.0]);
        let c = general_einsum(&a, &[0], &[], &b, &[0], &[]);
        assert!((c.get(&[0, 0]) - 10.0).abs() < 1e-10);
        assert!((c.get(&[1, 1]) - 21.0).abs() < 1e-10);
    }
    #[test]
    fn test_is_symmetric_modes_identity_matrix() {
        let mut t = DenseTensor::zeros(&[3, 3]);
        for i in 0..3 {
            t.set(&[i, i], 1.0);
        }
        assert!(is_symmetric_modes(&t, 0, 1, 1e-12));
    }
    #[test]
    fn test_is_symmetric_modes_antisymmetric() {
        let mut t = DenseTensor::zeros(&[2, 2]);
        t.set(&[0, 1], 1.0);
        t.set(&[1, 0], -1.0);
        assert!(!is_symmetric_modes(&t, 0, 1, 1e-12));
    }
    #[test]
    fn test_neo_hookean_stiffness_equals_isotropic() {
        let c1 = neo_hookean_stiffness(2.0, 3.0);
        let c2 = Tensor4::isotropic(2.0, 3.0);
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        assert!((c1.data[i][j][k][l] - c2.data[i][j][k][l]).abs() < 1e-12);
                    }
                }
            }
        }
    }
    #[test]
    fn test_cauchy_stress_linear_zero_strain() {
        let c = Tensor4::isotropic(1.0, 1.0);
        let eps = Tensor2::zero();
        let sigma = cauchy_stress_linear(&c, &eps);
        for i in 0..3 {
            for j in 0..3 {
                assert!(approx(sigma.data[i][j], 0.0));
            }
        }
    }
}
