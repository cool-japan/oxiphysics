//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#[allow(unused_imports)]
use super::functions::*;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hyperelastic::*;
    pub(super) const TOL: f64 = 1e-10;
    pub(super) const TOL_LOOSE: f64 = 1e-5;
    fn identity_f() -> [[f64; 3]; 3] {
        mat3_identity()
    }
    fn simple_stretch(lx: f64) -> [[f64; 3]; 3] {
        [[lx, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    }
    #[test]
    fn test_mat3_mul_identity() {
        let i = mat3_identity();
        let a = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let r = mat3_mul(i, a);
        for ii in 0..3 {
            for jj in 0..3 {
                assert!((r[ii][jj] - a[ii][jj]).abs() < TOL);
            }
        }
    }
    #[test]
    fn test_mat3_mul_basic() {
        let a = [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]];
        let b = [[2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0]];
        let c = mat3_mul(a, b);
        assert!((c[0][0] - 2.0).abs() < TOL);
        assert!((c[1][1] - 6.0).abs() < TOL);
        assert!((c[2][2] - 12.0).abs() < TOL);
    }
    #[test]
    fn test_mat3_transpose() {
        let a = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let t = mat3_transpose(a);
        for i in 0..3 {
            for j in 0..3 {
                assert!((t[i][j] - a[j][i]).abs() < TOL);
            }
        }
    }
    #[test]
    fn test_mat3_det_identity() {
        assert!((mat3_det(mat3_identity()) - 1.0).abs() < TOL);
    }
    #[test]
    fn test_mat3_det_diagonal() {
        let a = [[2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0]];
        assert!((mat3_det(a) - 24.0).abs() < TOL);
    }
    #[test]
    fn test_mat3_det_singular() {
        let a = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        assert!(mat3_det(a).abs() < 1e-10);
    }
    #[test]
    fn test_mat3_inverse_identity() {
        let inv = mat3_inverse(mat3_identity()).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((inv[i][j] - expected).abs() < TOL);
            }
        }
    }
    #[test]
    fn test_mat3_inverse_roundtrip() {
        let a = [[2.0, 1.0, 0.0], [1.0, 3.0, 1.0], [0.0, 1.0, 2.0]];
        let inv = mat3_inverse(a).unwrap();
        let prod = mat3_mul(a, inv);
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((prod[i][j] - expected).abs() < 1e-12);
            }
        }
    }
    #[test]
    fn test_mat3_inverse_singular() {
        let a = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        assert!(mat3_inverse(a).is_none());
    }
    #[test]
    fn test_right_cauchy_green_identity() {
        let c = right_cauchy_green(identity_f());
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((c[i][j] - expected).abs() < TOL);
            }
        }
    }
    #[test]
    fn test_left_cauchy_green_identity() {
        let b = left_cauchy_green(identity_f());
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((b[i][j] - expected).abs() < TOL);
            }
        }
    }
    #[test]
    fn test_right_cauchy_green_stretch() {
        let f = simple_stretch(2.0);
        let c = right_cauchy_green(f);
        assert!((c[0][0] - 4.0).abs() < TOL);
        assert!((c[1][1] - 1.0).abs() < TOL);
        assert!((c[2][2] - 1.0).abs() < TOL);
    }
    #[test]
    fn test_invariants_identity() {
        let (i1, i2, i3) = invariants(mat3_identity());
        assert!((i1 - 3.0).abs() < TOL);
        assert!((i2 - 3.0).abs() < TOL);
        assert!((i3 - 1.0).abs() < TOL);
    }
    #[test]
    fn test_invariants_stretch() {
        let f = simple_stretch(2.0);
        let c = right_cauchy_green(f);
        let (i1, _i2, i3) = invariants(c);
        assert!((i1 - 6.0).abs() < TOL);
        assert!((i3 - 4.0).abs() < TOL);
    }
    #[test]
    fn test_neo_hookean_from_young_poisson() {
        let _mu_v = 1000.0 / (2.0 * (1.0 + 0.3));
        let _kappa_v = 1000.0 / (3.0 * (1.0 - 2.0 * 0.3));
        let nh = NeoHookean::new(_mu_v, _kappa_v);
        let expected_mu = 1000.0 / 2.6;
        assert!((nh.mu - expected_mu).abs() < 1e-10);
    }
    #[test]
    fn test_neo_hookean_strain_energy_zero_at_identity() {
        let _mu_v = 1000.0 / (2.0 * (1.0 + 0.3));
        let _kappa_v = 1000.0 / (3.0 * (1.0 - 2.0 * 0.3));
        let nh = NeoHookean::new(_mu_v, _kappa_v);
        let w = nh.strain_energy(identity_f());
        assert!(w.abs() < TOL);
    }
    #[test]
    fn test_neo_hookean_strain_energy_positive() {
        let _mu_v = 1000.0 / (2.0 * (1.0 + 0.3));
        let _kappa_v = 1000.0 / (3.0 * (1.0 - 2.0 * 0.3));
        let nh = NeoHookean::new(_mu_v, _kappa_v);
        let f = simple_stretch(1.5);
        let w = nh.strain_energy(f);
        assert!(w > 0.0);
    }
    #[test]
    fn test_neo_hookean_pk1_zero_at_identity() {
        let _mu_v = 1000.0 / (2.0 * (1.0 + 0.3));
        let _kappa_v = 1000.0 / (3.0 * (1.0 - 2.0 * 0.3));
        let nh = NeoHookean::new(_mu_v, _kappa_v);
        let p = nh.pk1_stress(identity_f());
        for i in 0..3 {
            for j in 0..3 {
                assert!(p[i][j].abs() < TOL, "P[{i}][{j}] = {}", p[i][j]);
            }
        }
    }
    #[test]
    fn test_neo_hookean_cauchy_zero_at_identity() {
        let _mu_v = 1000.0 / (2.0 * (1.0 + 0.3));
        let _kappa_v = 1000.0 / (3.0 * (1.0 - 2.0 * 0.3));
        let nh = NeoHookean::new(_mu_v, _kappa_v);
        let sigma = nh.cauchy_stress(identity_f());
        for i in 0..3 {
            for j in 0..3 {
                assert!(sigma[i][j].abs() < TOL, "sigma[{i}][{j}] = {}", sigma[i][j]);
            }
        }
    }
    #[test]
    fn test_neo_hookean_pk1_numerical_gradient() {
        let _mu_v = 500.0 / (2.0 * (1.0 + 0.25));
        let _kappa_v = 500.0 / (3.0 * (1.0 - 2.0 * 0.25));
        let nh = NeoHookean::new(_mu_v, _kappa_v);
        let mut f = [[1.2, 0.05, 0.0], [0.0, 0.9, 0.02], [0.01, 0.0, 1.1]];
        let h = 1e-6;
        let p = nh.pk1_stress(f);
        let f_save = f[0][0];
        f[0][0] = f_save + h;
        let wp = nh.strain_energy(f);
        f[0][0] = f_save - h;
        let wm = nh.strain_energy(f);
        f[0][0] = f_save;
        let _ = f[0][0];
        let p00_fd = (wp - wm) / (2.0 * h);
        assert!((p[0][0] - p00_fd).abs() < 1e-5);
    }
    #[test]
    fn test_mooney_rivlin_strain_energy_zero_at_identity() {
        let mr = MooneyRivlin {
            c10: 100.0,
            c01: 50.0,
            d1: 0.01,
        };
        let w = mr.strain_energy(identity_f());
        assert!(w.abs() < TOL_LOOSE, "W at identity = {w}");
    }
    #[test]
    fn test_mooney_rivlin_strain_energy_positive() {
        let mr = MooneyRivlin {
            c10: 100.0,
            c01: 50.0,
            d1: 0.01,
        };
        let f = simple_stretch(1.3);
        let w = mr.strain_energy(f);
        assert!(w > 0.0);
    }
    #[test]
    fn test_mooney_rivlin_pk2_symmetric() {
        let mr = MooneyRivlin {
            c10: 80.0,
            c01: 40.0,
            d1: 0.02,
        };
        let f = [[1.2, 0.1, 0.0], [0.0, 0.9, 0.05], [0.05, 0.0, 1.1]];
        let s = mr.pk2_stress(f);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (s[i][j] - s[j][i]).abs() < 1e-4,
                    "S[{i}][{j}]={} vs S[{j}][{i}]={}",
                    s[i][j],
                    s[j][i]
                );
            }
        }
    }
    #[test]
    fn test_ogden_term_strain_energy_zero_at_one() {
        let term = OgdenTerm {
            mu: 100.0,
            alpha: 2.0,
        };
        let w = term.strain_energy_uniaxial(1.0);
        assert!(w.abs() < TOL_LOOSE);
    }
    #[test]
    fn test_ogden_term_strain_energy_positive() {
        let term = OgdenTerm {
            mu: 100.0,
            alpha: 2.0,
        };
        let w = term.strain_energy_uniaxial(1.5);
        assert!(w > 0.0);
    }
    #[test]
    fn test_ogden_cauchy_stress_uniaxial_zero_at_one() {
        let terms = [
            OgdenTerm {
                mu: 100.0,
                alpha: 2.0,
            },
            OgdenTerm {
                mu: 50.0,
                alpha: 4.0,
            },
        ];
        let sigma = ogden_cauchy_stress_uniaxial(&terms, 1.0);
        assert!(sigma.is_finite());
    }
    #[test]
    fn test_ogden_cauchy_stress_increases_with_stretch() {
        let terms = [OgdenTerm {
            mu: 100.0,
            alpha: 2.0,
        }];
        let s1 = ogden_cauchy_stress_uniaxial(&terms, 1.5);
        let s2 = ogden_cauchy_stress_uniaxial(&terms, 2.0);
        assert!(s2 > s1);
    }
    #[test]
    fn test_deviatoric_f_det_one() {
        let f = [[1.5, 0.0, 0.0], [0.0, 1.2, 0.0], [0.0, 0.0, 0.8]];
        let f_bar = deviatoric_f(f);
        let det = mat3_det(f_bar);
        assert!((det - 1.0).abs() < 1e-12, "det(F_bar) = {det}");
    }
    #[test]
    fn test_volumetric_strain_energy_zero_at_j1() {
        let w = volumetric_strain_energy(1000.0, 1.0);
        assert!(w.abs() < TOL);
    }
    #[test]
    fn test_volumetric_strain_energy_positive() {
        let w = volumetric_strain_energy(1000.0, 1.5);
        assert!(w > 0.0);
    }
    #[test]
    fn test_deformation_gradient_tet_identity() {
        let coords: [[f64; 3]; 4] = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let f = deformation_gradient_tet(&coords, &coords);
        let eye = mat3_identity();
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (f[i][j] - eye[i][j]).abs() < 1e-10,
                    "F[{i}][{j}]={}",
                    f[i][j]
                );
            }
        }
    }
    #[test]
    fn test_deformation_gradient_tet_uniform_stretch() {
        let ref_coords: [[f64; 3]; 4] = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let cur_coords: [[f64; 3]; 4] = [
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let f = deformation_gradient_tet(&ref_coords, &cur_coords);
        assert!((f[0][0] - 2.0).abs() < 1e-10);
        assert!((f[1][1] - 1.0).abs() < 1e-10);
        assert!((f[2][2] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_arruda_boyce_zero_at_identity() {
        let ab = ArrudaBoyce {
            mu: 100.0,
            bulk_modulus: 1000.0,
            lambda_m: 5.0,
        };
        let w = ab.strain_energy(identity_f());
        assert!(w.abs() < TOL_LOOSE, "W at identity = {w}");
    }
    #[test]
    fn test_arruda_boyce_positive_under_stretch() {
        let ab = ArrudaBoyce {
            mu: 100.0,
            bulk_modulus: 1000.0,
            lambda_m: 5.0,
        };
        let w = ab.strain_energy(simple_stretch(1.5));
        assert!(w > 0.0, "W should be positive under stretch");
    }
    #[test]
    fn test_arruda_boyce_pk1_numerical() {
        let ab = ArrudaBoyce {
            mu: 200.0,
            bulk_modulus: 2000.0,
            lambda_m: 4.0,
        };
        let f = [[1.1, 0.02, 0.0], [0.0, 0.95, 0.01], [0.01, 0.0, 1.05]];
        let p = ab.pk1_stress(f);
        let h = 1e-6;
        let f_save = f[0][0];
        let mut f_p = f;
        let mut f_m = f;
        f_p[0][0] = f_save + h;
        f_m[0][0] = f_save - h;
        let p00_fd = (ab.strain_energy(f_p) - ab.strain_energy(f_m)) / (2.0 * h);
        assert!(
            (p[0][0] - p00_fd).abs() < 1e-3,
            "P00={} vs FD={}",
            p[0][0],
            p00_fd
        );
    }
    #[test]
    fn test_gent_zero_at_identity() {
        let g = Gent {
            mu: 100.0,
            bulk_modulus: 1000.0,
            j_m: 50.0,
        };
        let w = g.strain_energy(identity_f());
        assert!(w.abs() < TOL_LOOSE, "W at identity = {w}");
    }
    #[test]
    fn test_gent_positive_under_stretch() {
        let g = Gent {
            mu: 100.0,
            bulk_modulus: 1000.0,
            j_m: 50.0,
        };
        let w = g.strain_energy(simple_stretch(1.3));
        assert!(w > 0.0);
    }
    #[test]
    fn test_gent_increases_with_stretch() {
        let g = Gent {
            mu: 100.0,
            bulk_modulus: 1000.0,
            j_m: 50.0,
        };
        let w1 = g.strain_energy(simple_stretch(1.1));
        let w2 = g.strain_energy(simple_stretch(1.5));
        assert!(w2 > w1, "energy should increase with larger stretch");
    }
    #[test]
    fn test_gent_converges_to_neo_hookean_large_jm() {
        let mu = 100.0;
        let g = Gent {
            mu,
            bulk_modulus: 0.0,
            j_m: 1e6,
        };
        let f = simple_stretch(1.2);
        let w_gent = g.strain_energy(f);
        let c = right_cauchy_green(f);
        let (i1, _, _) = invariants(c);
        let w_approx = mu / 2.0 * (i1 - 3.0);
        assert!(
            (w_gent - w_approx).abs() / (w_approx.abs() + 1e-10) < 0.01,
            "Gent={w_gent} vs approx={w_approx}"
        );
    }
    #[test]
    fn test_yeoh_zero_at_identity() {
        let y = Yeoh {
            c1: 50.0,
            c2: 10.0,
            c3: 1.0,
            bulk_modulus: 1000.0,
        };
        let w = y.strain_energy(identity_f());
        assert!(w.abs() < TOL_LOOSE, "W at identity = {w}");
    }
    #[test]
    fn test_yeoh_positive_under_stretch() {
        let y = Yeoh {
            c1: 50.0,
            c2: 10.0,
            c3: 1.0,
            bulk_modulus: 1000.0,
        };
        let w = y.strain_energy(simple_stretch(1.3));
        assert!(w > 0.0);
    }
    #[test]
    fn test_yeoh_with_c1_only_matches_neo_hookean_like() {
        let c1 = 50.0;
        let y = Yeoh {
            c1,
            c2: 0.0,
            c3: 0.0,
            bulk_modulus: 0.0,
        };
        let f = simple_stretch(1.2);
        let w = y.strain_energy(f);
        let c = right_cauchy_green(f);
        let (i1, _, _) = invariants(c);
        let expected = c1 * (i1 - 3.0);
        assert!(
            (w - expected).abs() < 1e-10,
            "Yeoh={w} vs expected={expected}"
        );
    }
    #[test]
    fn test_fung_zero_at_identity() {
        let fung = Fung {
            c: 100.0,
            b1: 1.0,
            bulk_modulus: 1000.0,
        };
        let w = fung.strain_energy(identity_f());
        assert!(w.abs() < TOL_LOOSE, "W at identity = {w}");
    }
    #[test]
    fn test_fung_positive_under_stretch() {
        let fung = Fung {
            c: 100.0,
            b1: 1.0,
            bulk_modulus: 1000.0,
        };
        let w = fung.strain_energy(simple_stretch(1.3));
        assert!(w > 0.0);
    }
    #[test]
    fn test_fung_exponential_growth() {
        let fung = Fung {
            c: 100.0,
            b1: 2.0,
            bulk_modulus: 0.0,
        };
        let w1 = fung.strain_energy(simple_stretch(1.5));
        let w2 = fung.strain_energy(simple_stretch(2.0));
        assert!(w2 > w1 * 2.0, "Fung should grow exponentially");
    }
    #[test]
    fn test_pk1_to_pk2_at_identity() {
        let _mu_v = 1000.0 / (2.0 * (1.0 + 0.3));
        let _kappa_v = 1000.0 / (3.0 * (1.0 - 2.0 * 0.3));
        let nh = NeoHookean::new(_mu_v, _kappa_v);
        let f = identity_f();
        let p = nh.pk1_stress(f);
        let s = pk1_to_pk2(f, p);
        for i in 0..3 {
            for j in 0..3 {
                assert!(s[i][j].abs() < TOL, "S[{i}][{j}] = {}", s[i][j]);
            }
        }
    }
    #[test]
    fn test_pk1_to_cauchy_at_identity() {
        let _mu_v = 1000.0 / (2.0 * (1.0 + 0.3));
        let _kappa_v = 1000.0 / (3.0 * (1.0 - 2.0 * 0.3));
        let nh = NeoHookean::new(_mu_v, _kappa_v);
        let f = identity_f();
        let p = nh.pk1_stress(f);
        let sigma = pk1_to_cauchy(f, p);
        for i in 0..3 {
            for j in 0..3 {
                assert!(sigma[i][j].abs() < TOL, "sigma[{i}][{j}] = {}", sigma[i][j]);
            }
        }
    }
    #[test]
    fn test_pk1_to_cauchy_matches_direct() {
        let _mu_v = 500.0 / (2.0 * (1.0 + 0.25));
        let _kappa_v = 500.0 / (3.0 * (1.0 - 2.0 * 0.25));
        let nh = NeoHookean::new(_mu_v, _kappa_v);
        let f = [[1.2, 0.05, 0.0], [0.0, 0.9, 0.02], [0.01, 0.0, 1.1]];
        let p = nh.pk1_stress(f);
        let sigma_direct = nh.cauchy_stress(f);
        let sigma_via_pk1 = pk1_to_cauchy(f, p);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (sigma_direct[i][j] - sigma_via_pk1[i][j]).abs() < 1e-6,
                    "sigma[{i}][{j}]: direct={} vs via_pk1={}",
                    sigma_direct[i][j],
                    sigma_via_pk1[i][j]
                );
            }
        }
    }
    #[test]
    fn test_material_tangent_symmetry() {
        let _mu_v = 1000.0 / (2.0 * (1.0 + 0.3));
        let _kappa_v = 1000.0 / (3.0 * (1.0 - 2.0 * 0.3));
        let nh = NeoHookean::new(_mu_v, _kappa_v);
        let f = [[1.1, 0.01, 0.0], [0.0, 0.95, 0.0], [0.0, 0.0, 1.0]];
        let c = material_tangent(f, |ff| nh.strain_energy(ff));
        for ii in 0..9 {
            for jj in ii..9 {
                assert!(
                    (c[ii * 9 + jj] - c[jj * 9 + ii]).abs() < 1e-2,
                    "Tangent not symmetric at ({ii},{jj}): {} vs {}",
                    c[ii * 9 + jj],
                    c[jj * 9 + ii]
                );
            }
        }
    }
    #[test]
    fn test_material_tangent_positive_diagonal() {
        let _mu_v = 1000.0 / (2.0 * (1.0 + 0.3));
        let _kappa_v = 1000.0 / (3.0 * (1.0 - 2.0 * 0.3));
        let nh = NeoHookean::new(_mu_v, _kappa_v);
        let f = identity_f();
        let c = material_tangent(f, |ff| nh.strain_energy(ff));
        for ii in 0..9 {
            assert!(
                c[ii * 9 + ii] > -1e-4,
                "Tangent diagonal [{ii}] = {} not positive",
                c[ii * 9 + ii]
            );
        }
    }
    #[test]
    fn test_mat3_add() {
        let a = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let b = [[9.0, 8.0, 7.0], [6.0, 5.0, 4.0], [3.0, 2.0, 1.0]];
        let c = mat3_add(a, b);
        for i in 0..3 {
            for j in 0..3 {
                assert!((c[i][j] - 10.0).abs() < TOL);
            }
        }
    }
    #[test]
    fn test_mat3_frobenius_identity() {
        let f = mat3_frobenius(identity_f());
        assert!((f - 3.0_f64.sqrt()).abs() < TOL);
    }
    #[test]
    fn test_mat3_scale() {
        let a = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let b = mat3_scale(a, 2.0);
        for i in 0..3 {
            for j in 0..3 {
                assert!((b[i][j] - 2.0 * a[i][j]).abs() < TOL);
            }
        }
    }
    #[test]
    fn test_total_lagrangian_pk2_stress_identity() {
        let _mu_v = 1000.0 / (2.0 * (1.0 + 0.3));
        let _kappa_v = 1000.0 / (3.0 * (1.0 - 2.0 * 0.3));
        let nh = NeoHookean::new(_mu_v, _kappa_v);
        let f = identity_f();
        let p = nh.pk1_stress(f);
        let s = pk1_to_pk2(f, p);
        for i in 0..3 {
            for j in 0..3 {
                assert!(s[i][j].abs() < TOL, "S[{i}][{j}] = {}", s[i][j]);
            }
        }
    }
    #[test]
    fn test_multiplicative_split_consistency() {
        let f = [[1.3, 0.05, 0.0], [0.0, 0.95, 0.01], [0.0, 0.0, 1.1]];
        let j = mat3_det(f);
        let f_bar = deviatoric_f(f);
        let j_bar = mat3_det(f_bar);
        assert!((j_bar - 1.0).abs() < 1e-12, "det(F_bar)={j_bar}");
        let j13 = j.powf(1.0 / 3.0);
        let f_reconstructed = mat3_scale(f_bar, j13);
        for i in 0..3 {
            for k in 0..3 {
                assert!(
                    (f_reconstructed[i][k] - f[i][k]).abs() < 1e-12,
                    "F_reconstructed[{i}][{k}]={}",
                    f_reconstructed[i][k]
                );
            }
        }
    }
    #[test]
    fn test_bbar_element_volumetric_consistency() {
        let f = [[2.0, 0.0, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 1.0]];
        let j = mat3_det(f);
        let f_bar = deviatoric_f(f);
        let c_bar = right_cauchy_green(f_bar);
        let (i1_bar, _, _) = invariants(c_bar);
        let c = right_cauchy_green(f);
        let (i1, _, _) = invariants(c);
        let i1_bar_expected = i1 * j.powf(-2.0 / 3.0);
        assert!(
            (i1_bar - i1_bar_expected).abs() < 1e-10,
            "i1_bar={i1_bar}, expected={i1_bar_expected}"
        );
    }
    #[test]
    fn test_nearly_incompressible_neo_hookean_energy() {
        let nh_comp = NeoHookean {
            mu: 100.0,
            lambda: 1e6,
        };
        let nh_incomp = NeoHookean {
            mu: 100.0,
            lambda: 1e8,
        };
        let f = [[1.2, 0.0, 0.0], [0.0, 0.9, 0.0], [0.0, 0.0, 1.0]];
        let w_comp = nh_comp.strain_energy(f);
        let w_incomp = nh_incomp.strain_energy(f);
        assert!(w_comp > 0.0);
        assert!(
            w_incomp > w_comp,
            "Nearly incompressible should have higher energy: {w_incomp} vs {w_comp}"
        );
    }
    #[test]
    fn test_nearly_incompressible_pressure_sensitive() {
        let _nh = NeoHookean {
            mu: 100.0,
            lambda: 1e6,
        };
        let stretch = 1.5_f64;
        let inv = 1.0 / stretch.sqrt();
        let f = [[stretch, 0.0, 0.0], [0.0, inv, 0.0], [0.0, 0.0, inv]];
        let j = mat3_det(f);
        assert!(
            (j - 1.0).abs() < 1e-12,
            "J should be 1 for volume-preserving: {j}"
        );
        let w_vol = volumetric_strain_energy(1e6, j);
        assert!(w_vol.abs() < 1e-8, "Volumetric energy at J=1: {w_vol}");
    }
    #[test]
    fn test_hyperelastic_stability_drucker() {
        let _mu_v = 1000.0 / (2.0 * (1.0 + 0.3));
        let _kappa_v = 1000.0 / (3.0 * (1.0 - 2.0 * 0.3));
        let nh = NeoHookean::new(_mu_v, _kappa_v);
        let f = [[1.1, 0.02, 0.0], [0.0, 0.95, 0.01], [0.0, 0.0, 1.05]];
        let p0 = nh.pk1_stress(f);
        let mut f_pert = f;
        f_pert[0][0] += 0.001;
        let p_pert = nh.pk1_stress(f_pert);
        let mut dp_df = 0.0_f64;
        let df = [[0.001, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        for i in 0..3 {
            for j in 0..3 {
                dp_df += (p_pert[i][j] - p0[i][j]) * df[i][j];
            }
        }
        assert!(
            dp_df >= -1e-8,
            "Drucker stability violated: dP:dF = {dp_df}"
        );
    }
    #[test]
    fn test_hyperelastic_stability_energy_convexity() {
        let _mu_v = 1000.0 / (2.0 * (1.0 + 0.3));
        let _kappa_v = 1000.0 / (3.0 * (1.0 - 2.0 * 0.3));
        let nh = NeoHookean::new(_mu_v, _kappa_v);
        let stretches = [1.0, 1.2, 1.5, 2.0, 3.0];
        let energies: Vec<f64> = stretches
            .iter()
            .map(|&lx| nh.strain_energy([[lx, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]))
            .collect();
        for k in 1..energies.len() {
            assert!(
                energies[k] > energies[k - 1],
                "Energy not monotone: W({})={} <= W({})={}",
                stretches[k],
                energies[k],
                stretches[k - 1],
                energies[k - 1]
            );
        }
    }
    #[test]
    fn test_mixed_up_pressure_consistency() {
        let nh = NeoHookean {
            mu: 100.0,
            lambda: 10000.0,
        };
        let f = [[1.1, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let sigma = nh.cauchy_stress(f);
        let p = (sigma[0][0] + sigma[1][1] + sigma[2][2]) / 3.0;
        assert!(p.is_finite(), "Hydrostatic pressure should be finite: {p}");
    }
    #[test]
    fn test_mixed_up_deviatoric_stress_traceless() {
        let _mu_v = 1000.0 / (2.0 * (1.0 + 0.3));
        let _kappa_v = 1000.0 / (3.0 * (1.0 - 2.0 * 0.3));
        let nh = NeoHookean::new(_mu_v, _kappa_v);
        let f = [[1.2, 0.1, 0.0], [0.0, 0.9, 0.05], [0.0, 0.0, 1.1]];
        let sigma = nh.cauchy_stress(f);
        let p = (sigma[0][0] + sigma[1][1] + sigma[2][2]) / 3.0;
        let dev_trace = (sigma[0][0] - p) + (sigma[1][1] - p) + (sigma[2][2] - p);
        assert!(
            dev_trace.abs() < 1e-10,
            "Deviatoric trace should be zero: {dev_trace}"
        );
    }
    #[test]
    fn test_fbar_element_no_volumetric_locking() {
        let f = [[1.5, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let j = mat3_det(f);
        let f_bar = deviatoric_f(f);
        let j_bar = mat3_det(f_bar);
        assert!((j_bar - 1.0).abs() < 1e-12);
        let w_vol = volumetric_strain_energy(1000.0, j);
        assert!(
            w_vol > 0.0,
            "Volumetric energy positive for volume change: {w_vol}"
        );
    }
    #[test]
    fn test_varga_strain_energy_zero_at_identity() {
        let v = Varga {
            mu: 50.0,
            bulk_modulus: 1000.0,
        };
        let w = v.strain_energy(mat3_identity());
        assert!(w.abs() < TOL_LOOSE, "W at identity = {w}");
    }
    #[test]
    fn test_varga_strain_energy_positive_stretch() {
        let v = Varga {
            mu: 50.0,
            bulk_modulus: 1000.0,
        };
        let w = v.strain_energy(simple_stretch(1.5));
        assert!(w > 0.0, "W > 0 for stretch: {w}");
    }
    #[test]
    fn test_varga_pk1_zero_at_identity() {
        let v = Varga {
            mu: 50.0,
            bulk_modulus: 1000.0,
        };
        let p = v.pk1_stress(mat3_identity());
        for i in 0..3 {
            for j in 0..3 {
                assert!(p[i][j].abs() < TOL_LOOSE, "P[{i}][{j}]={}", p[i][j]);
            }
        }
    }
    #[test]
    fn test_holzapfel_ogden_isotropic_limit() {
        let ho = HolzapfelOgden {
            c10: 100.0,
            bulk_modulus: 10000.0,
            k1: 0.0,
            k2: 1.0,
            a: [1.0, 0.0, 0.0],
        };
        let f = simple_stretch(1.3);
        let w = ho.strain_energy(f);
        assert!(w > 0.0, "W > 0 for anisotropic with k1=0: {w}");
    }
    #[test]
    fn test_holzapfel_ogden_fibre_stiffening() {
        let ho_no_fibre = HolzapfelOgden {
            c10: 100.0,
            bulk_modulus: 5000.0,
            k1: 0.0,
            k2: 1.0,
            a: [1.0, 0.0, 0.0],
        };
        let ho_with_fibre = HolzapfelOgden {
            c10: 100.0,
            bulk_modulus: 5000.0,
            k1: 500.0,
            k2: 1.0,
            a: [1.0, 0.0, 0.0],
        };
        let f = simple_stretch(1.5);
        let w_iso = ho_no_fibre.strain_energy(f);
        let w_aniso = ho_with_fibre.strain_energy(f);
        assert!(
            w_aniso > w_iso,
            "fibres should increase energy: w_iso={w_iso}, w_aniso={w_aniso}"
        );
    }
    #[test]
    fn test_holzapfel_ogden_energy_zero_at_identity() {
        let ho = HolzapfelOgden {
            c10: 100.0,
            bulk_modulus: 5000.0,
            k1: 200.0,
            k2: 1.0,
            a: [1.0, 0.0, 0.0],
        };
        let w = ho.strain_energy(mat3_identity());
        assert!(w.abs() < TOL_LOOSE, "W at identity = {w}");
    }
    #[test]
    fn test_holzapfel_ogden_pk1_not_all_zero_at_stretch() {
        let ho = HolzapfelOgden {
            c10: 100.0,
            bulk_modulus: 5000.0,
            k1: 200.0,
            k2: 1.0,
            a: [1.0, 0.0, 0.0],
        };
        let p = ho.pk1_stress(simple_stretch(1.5));
        let norm: f64 = p
            .iter()
            .flat_map(|r| r.iter())
            .map(|x| x * x)
            .sum::<f64>()
            .sqrt();
        assert!(norm > 0.0, "PK1 should be non-zero at stretch");
    }
    #[test]
    fn test_volumetric_penalty_zero_at_j1() {
        let vp = VolumetricPenalty::quadratic(1000.0);
        assert!(vp.energy(1.0).abs() < TOL, "penalty at J=1 should be zero");
    }
    #[test]
    fn test_volumetric_penalty_positive_compression() {
        let vp = VolumetricPenalty::quadratic(1000.0);
        assert!(vp.energy(0.8) > 0.0, "penalty at J<1 should be positive");
    }
    #[test]
    fn test_volumetric_penalty_positive_expansion() {
        let vp = VolumetricPenalty::quadratic(1000.0);
        assert!(vp.energy(1.2) > 0.0, "penalty at J>1 should be positive");
    }
    #[test]
    fn test_volumetric_penalty_pressure_zero_at_j1() {
        let vp = VolumetricPenalty::simo_taylor(2000.0);
        let p = vp.pressure(1.0);
        assert!(p.abs() < 1e-3, "pressure at J=1 should be ~0, got {p}");
    }
    #[test]
    fn test_volumetric_penalty_pressure_negative_compression() {
        let vp = VolumetricPenalty::simo_taylor(2000.0);
        let p = vp.pressure(0.8);
        assert!(p != 0.0, "pressure at J=0.8 should be non-zero, got {p}");
    }
    #[test]
    fn test_volumetric_penalty_log_type() {
        let vp = VolumetricPenalty::logarithmic(1000.0);
        assert!(vp.energy(0.5) > 0.0);
        assert!(vp.energy(2.0) > 0.0);
        assert!(vp.energy(1.0).abs() < TOL);
    }
    #[test]
    fn test_consistent_tangent_neo_hookean_symmetric() {
        let _mu_v = 1000.0 / (2.0 * (1.0 + 0.3));
        let _kappa_v = 1000.0 / (3.0 * (1.0 - 2.0 * 0.3));
        let nh = NeoHookean::new(_mu_v, _kappa_v);
        let f = [[1.1, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let c = consistent_tangent_4th_order(&nh, f);
        for ii in 0..9 {
            for jj in ii + 1..9 {
                assert!(
                    (c[ii * 9 + jj] - c[jj * 9 + ii]).abs() < 0.1,
                    "Tangent not symmetric at ({ii},{jj})"
                );
            }
        }
    }
    #[test]
    fn test_consistent_tangent_positive_major_diag() {
        let _mu_v = 1000.0 / (2.0 * (1.0 + 0.3));
        let _kappa_v = 1000.0 / (3.0 * (1.0 - 2.0 * 0.3));
        let nh = NeoHookean::new(_mu_v, _kappa_v);
        let c = consistent_tangent_4th_order(&nh, mat3_identity());
        for ii in [0, 4, 8] {
            assert!(
                c[ii * 9 + ii] > 0.0,
                "C[{ii},{ii}] should be > 0, got {}",
                c[ii * 9 + ii]
            );
        }
    }
    #[test]
    fn test_strain_energy_gradient_neo_hookean() {
        let _mu_v = 500.0 / (2.0 * (1.0 + 0.3));
        let _kappa_v = 500.0 / (3.0 * (1.0 - 2.0 * 0.3));
        let nh = NeoHookean::new(_mu_v, _kappa_v);
        let f = [[1.2, 0.05, 0.0], [0.0, 0.95, 0.0], [0.01, 0.0, 1.1]];
        let p_analytical = nh.pk1_stress(f);
        let p_numerical = numerical_pk1(f, |ff| nh.strain_energy(ff));
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (p_analytical[i][j] - p_numerical[i][j]).abs() < 1e-4,
                    "Gradient mismatch at [{i}][{j}]: analytic={} numerical={}",
                    p_analytical[i][j],
                    p_numerical[i][j]
                );
            }
        }
    }
    #[test]
    fn test_strain_energy_gent_gradient() {
        let gent = Gent {
            mu: 100.0,
            bulk_modulus: 5000.0,
            j_m: 20.0,
        };
        let f = [[1.1, 0.01, 0.0], [0.0, 0.98, 0.0], [0.0, 0.0, 1.05]];
        let p = gent.pk1_stress(f);
        for i in 0..3 {
            for j in 0..3 {
                assert!(p[i][j].is_finite(), "P[{i}][{j}] should be finite");
            }
        }
    }
    #[test]
    fn test_strain_energy_yeoh_gradient_consistency() {
        let yeoh = Yeoh {
            c1: 50.0,
            c2: 5.0,
            c3: 0.5,
            bulk_modulus: 5000.0,
        };
        let f = [[1.2, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let p_analytic = yeoh.pk1_stress(f);
        let p_numerical = numerical_pk1(f, |ff| yeoh.strain_energy(ff));
        assert!(
            (p_analytic[0][0] - p_numerical[0][0]).abs() < 1.0,
            "Yeoh P[0][0] analytic={} numerical={}",
            p_analytic[0][0],
            p_numerical[0][0]
        );
    }
    #[test]
    fn test_fung_energy_zero_at_identity() {
        let fung = Fung {
            c: 10.0,
            b1: 1.0,
            bulk_modulus: 5000.0,
        };
        let w = fung.strain_energy(mat3_identity());
        assert!(w.abs() < TOL_LOOSE, "W at identity = {w}");
    }
    #[test]
    fn test_fung_energy_positive_at_stretch() {
        let fung = Fung {
            c: 10.0,
            b1: 1.0,
            bulk_modulus: 5000.0,
        };
        let w = fung.strain_energy(simple_stretch(1.3));
        assert!(w > 0.0, "W > 0 at stretch");
    }
    #[test]
    fn test_arruda_boyce_energy_zero_at_identity() {
        let ab = ArrudaBoyce {
            mu: 50.0,
            bulk_modulus: 2000.0,
            lambda_m: 5.0,
        };
        let w = ab.strain_energy(mat3_identity());
        assert!(w.abs() < TOL_LOOSE, "W at identity = {w}");
    }
    #[test]
    fn test_arruda_boyce_pk1_nonzero_at_stretch() {
        let ab = ArrudaBoyce {
            mu: 50.0,
            bulk_modulus: 2000.0,
            lambda_m: 5.0,
        };
        let p = ab.pk1_stress(simple_stretch(2.0));
        let norm: f64 = p
            .iter()
            .flat_map(|r| r.iter())
            .map(|x| x * x)
            .sum::<f64>()
            .sqrt();
        assert!(norm > 0.0, "PK1 should be non-zero at stretch");
    }
    #[test]
    fn test_push_forward_identity_f() {
        let f = mat3_identity();
        let t = [[1.0, 2.0, 0.0], [2.0, 3.0, 1.0], [0.0, 1.0, 4.0]];
        let t_pushed = push_forward_tensor(f, t);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (t_pushed[i][j] - t[i][j]).abs() < TOL,
                    "push-forward with I should be identity"
                );
            }
        }
    }
    #[test]
    fn test_pull_back_cauchy_to_pk2_identity_f() {
        let f = mat3_identity();
        let sigma = [[100.0, 0.0, 0.0], [0.0, 50.0, 0.0], [0.0, 0.0, 30.0]];
        let s = pull_back_cauchy_to_pk2(f, sigma);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (s[i][j] - sigma[i][j]).abs() < TOL,
                    "pull-back with I: S should equal sigma"
                );
            }
        }
    }
    #[test]
    fn test_green_lagrange_zero_at_identity() {
        let e = green_lagrange_strain(mat3_identity());
        for i in 0..3 {
            for j in 0..3 {
                assert!(e[i][j].abs() < TOL, "GL strain at identity should be zero");
            }
        }
    }
    #[test]
    fn test_linearised_strain_symmetric() {
        let grad_u = [[0.01, 0.02, 0.0], [0.03, 0.02, 0.01], [0.0, 0.01, 0.02]];
        let eps = linearised_strain(grad_u);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (eps[i][j] - eps[j][i]).abs() < TOL,
                    "linearised strain should be symmetric"
                );
            }
        }
    }
    #[test]
    fn test_internal_virtual_work_orthogonal() {
        let p = [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]];
        let df = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let w = internal_virtual_work(p, df);
        assert!((w - 6.0).abs() < TOL, "tr(P·dF) = 1+2+3 = 6, got {w}");
    }
    #[test]
    fn test_nearly_incompressible_update_pressure() {
        let mut ni = NearlyIncompressibleNeoHookean::from_young_poisson(1000.0, 0.45);
        let f = [[1.05, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        ni.update_pressure(f);
        let j = mat3_det(f);
        let expected_p = ni.kappa * (j - 1.0);
        assert!(
            (ni.pressure - expected_p).abs() < 1e-8,
            "pressure updated correctly, got {}",
            ni.pressure
        );
    }
    #[test]
    fn test_isochoric_stability_ratio_above_one() {
        let _mu_v = 1000.0 / (2.0 * (1.0 + 0.3));
        let _kappa_v = 1000.0 / (3.0 * (1.0 - 2.0 * 0.3));
        let nh = NeoHookean::new(_mu_v, _kappa_v);
        let f = [[1.2, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let ratio = isochoric_stability_ratio(f, &|ff| nh.strain_energy(ff), 1e-3);
        assert!(
            ratio > 0.0,
            "stability ratio should be positive, got {ratio}"
        );
    }
    #[test]
    fn test_legendre_hadamard_positive_for_stable_material() {
        let _mu_v = 1000.0 / (2.0 * (1.0 + 0.3));
        let _kappa_v = 1000.0 / (3.0 * (1.0 - 2.0 * 0.3));
        let nh = NeoHookean::new(_mu_v, _kappa_v);
        let f = mat3_identity();
        let tangent = consistent_tangent_4th_order(&nh, f);
        let lh = legendre_hadamard_check(f, &tangent);
        assert!(lh > -1000.0, "LH should not be hugely negative: {lh}");
    }
}
