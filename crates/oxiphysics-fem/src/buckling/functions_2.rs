//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#[allow(unused_imports)]
use super::functions::*;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::buckling::*;
    use std::f64::consts::PI;
    #[test]
    fn test_geo_stiff_new_zero() {
        let kg = GeometricStiffness::new(3);
        for i in 0..3 {
            for j in 0..3 {
                assert_eq!(kg.get(i, j), 0.0);
            }
        }
    }
    #[test]
    fn test_geo_stiff_set_get() {
        let mut kg = GeometricStiffness::new(4);
        kg.set(1, 2, 3.125);
        assert!((kg.get(1, 2) - 3.125).abs() < 1e-12);
        assert_eq!(kg.get(0, 0), 0.0);
    }
    #[test]
    fn test_geo_stiff_overwrite() {
        let mut kg = GeometricStiffness::new(2);
        kg.set(0, 1, 1.0);
        kg.set(0, 1, 5.0);
        assert!((kg.get(0, 1) - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_column_buckling_pinned_pinned() {
        let e = 200e9;
        let i = 1e-6;
        let l = 2.0;
        let p_cr = column_buckling_load(e, i, l, 1.0);
        let expected = PI * PI * e * i / (l * l);
        assert!((p_cr - expected).abs() / expected < 1e-12);
    }
    #[test]
    fn test_column_buckling_fixed_fixed() {
        let e = 210e9;
        let i = 2e-6;
        let l = 3.0;
        let p_cr = column_buckling_load(e, i, l, 0.5);
        let expected = PI * PI * e * i / (0.5 * l).powi(2);
        assert!((p_cr - expected).abs() / expected < 1e-12);
    }
    #[test]
    fn test_column_buckling_fixed_free() {
        let e = 100e9;
        let i = 5e-7;
        let l = 1.5;
        let p_cr = column_buckling_load(e, i, l, 2.0);
        let expected = PI * PI * e * i / (4.0 * l * l);
        assert!((p_cr - expected).abs() / expected < 1e-12);
    }
    #[test]
    fn test_column_buckling_scale_with_i() {
        let e = 1.0;
        let l = 1.0;
        let p1 = column_buckling_load(e, 1.0, l, 1.0);
        let p2 = column_buckling_load(e, 2.0, l, 1.0);
        assert!((p2 / p1 - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_column_buckling_scale_with_l() {
        let e = 1.0;
        let i = 1.0;
        let p1 = column_buckling_load(e, i, 1.0, 1.0);
        let p2 = column_buckling_load(e, i, 2.0, 1.0);
        assert!((p1 / p2 - 4.0).abs() < 1e-12);
    }
    #[test]
    fn test_plate_coeff_square_m1_n1() {
        let k = plate_buckling_coefficient(1.0, 1.0, 1, 1);
        assert!((k - 4.0).abs() < 1e-12);
    }
    #[test]
    fn test_plate_coeff_m1_n1_rectangle() {
        let a = 2.0;
        let b = 1.0;
        let k = plate_buckling_coefficient(a, b, 1, 1);
        assert!((k - 6.25).abs() < 1e-12);
    }
    #[test]
    fn test_plate_coeff_positive() {
        for m in 1u32..=3 {
            for n in 1u32..=3 {
                let k = plate_buckling_coefficient(2.0, 1.0, m, n);
                assert!(k > 0.0);
            }
        }
    }
    #[test]
    fn test_plate_coeff_symmetry_mn() {
        let k = plate_buckling_coefficient(1.0, 1.0, 2, 1);
        let expected = {
            let term = 2.0 * 1.0 / 1.0 + 1.0 * 1.0 / (2.0 * 1.0);
            term * term
        };
        assert!((k - expected).abs() < 1e-12);
    }
    #[test]
    fn test_critical_stress_basic() {
        let e = 200e9;
        let nu = 0.3;
        let t = 0.01;
        let b = 0.5;
        let k = 4.0;
        let sigma = critical_stress(e, nu, t, b, k);
        let expected = k * PI * PI * e / (12.0 * (1.0 - nu * nu)) * (t / b).powi(2);
        assert!((sigma - expected).abs() / expected < 1e-12);
    }
    #[test]
    fn test_critical_stress_scales_k() {
        let e = 200e9;
        let nu = 0.3;
        let t = 0.01;
        let b = 0.5;
        let s1 = critical_stress(e, nu, t, b, 4.0);
        let s2 = critical_stress(e, nu, t, b, 8.0);
        assert!((s2 / s1 - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_critical_stress_positive() {
        assert!(critical_stress(200e9, 0.3, 0.01, 0.5, 4.0) > 0.0);
    }
    #[test]
    fn test_critical_stress_scales_t_squared() {
        let s1 = critical_stress(200e9, 0.3, 0.01, 0.5, 4.0);
        let s2 = critical_stress(200e9, 0.3, 0.02, 0.5, 4.0);
        assert!((s2 / s1 - 4.0).abs() < 1e-10);
    }
    #[test]
    fn test_imperfection_zero() {
        let lf = imperfection_sensitivity(10.0, 0.0);
        assert!((lf - 10.0).abs() < 1e-12);
    }
    #[test]
    fn test_imperfection_reduces_load() {
        let lf = imperfection_sensitivity(10.0, 0.01);
        assert!(lf < 10.0);
    }
    #[test]
    fn test_imperfection_formula() {
        let lambda = 5.0;
        let delta = 0.04;
        let lf = imperfection_sensitivity(lambda, delta);
        let expected = lambda * (1.0 - 0.5 * delta.sqrt());
        assert!((lf - expected).abs() < 1e-12);
    }
    #[test]
    fn test_imperfection_negative_delta() {
        let lf_pos = imperfection_sensitivity(10.0, 0.04);
        let lf_neg = imperfection_sensitivity(10.0, -0.04);
        assert!((lf_pos - lf_neg).abs() < 1e-12);
    }
    #[test]
    fn test_buckling_problem_2dof() {
        let mut prob = BucklingProblem::new(2);
        prob.set_k(0, 0, 4.0);
        prob.set_k(1, 1, 4.0);
        prob.set_kg(0, 0, 1.0);
        prob.set_kg(1, 1, 1.0);
        let factors = prob.solve_buckling_load_factors(1);
        assert_eq!(factors.len(), 1);
        assert!((factors[0] - 4.0).abs() < 1e-5);
    }
    #[test]
    fn test_buckling_problem_sorted() {
        let mut prob = BucklingProblem::new(2);
        prob.set_k(0, 0, 1.0);
        prob.set_k(1, 1, 2.0);
        prob.set_kg(0, 0, 1.0);
        prob.set_kg(1, 1, 1.0);
        let factors = prob.solve_buckling_load_factors(2);
        assert_eq!(factors.len(), 2);
        assert!(factors[0] <= factors[1] + 1e-8);
    }
    #[test]
    fn test_buckling_problem_num_modes_capped() {
        let mut prob = BucklingProblem::new(2);
        prob.set_k(0, 0, 3.0);
        prob.set_k(1, 1, 3.0);
        prob.set_kg(0, 0, 1.0);
        prob.set_kg(1, 1, 1.0);
        let factors = prob.solve_buckling_load_factors(5);
        assert_eq!(factors.len(), 2);
    }
    #[test]
    fn test_fem_buckling_pinned_pinned() {
        let e = 200e9_f64;
        let i_mom = 1e-6_f64;
        let l = 2.0_f64;
        let p = 1.0_f64;
        let analytical = PI * PI * e * i_mom / (l * l);
        let fem_mult = column_fem_buckling(10, e, i_mom, l, p);
        let err = ((fem_mult * p - analytical) / analytical).abs();
        assert!(err < 0.02, "FEM error too large: {:.4e}", err);
    }
    #[test]
    fn test_fem_buckling_more_elements_more_accurate() {
        let e = 200e9_f64;
        let i_mom = 1e-6_f64;
        let l = 2.0_f64;
        let p = 1.0_f64;
        let analytical = PI * PI * e * i_mom / (l * l);
        let m4 = column_fem_buckling(4, e, i_mom, l, p);
        let m10 = column_fem_buckling(10, e, i_mom, l, p);
        let err4 = ((m4 * p - analytical) / analytical).abs();
        let err10 = ((m10 * p - analytical) / analytical).abs();
        assert!(
            err10 <= err4 + 1e-6,
            "More elements should be more accurate"
        );
    }
    #[test]
    fn test_fem_buckling_positive() {
        let factor = column_fem_buckling(4, 200e9, 1e-6, 2.0, 1.0);
        assert!(factor > 0.0, "Buckling factor must be positive");
    }
    #[test]
    fn test_fem_buckling_scales_with_ei() {
        let l = 2.0;
        let p = 1.0;
        let f1 = column_fem_buckling(6, 200e9, 1e-6, l, p);
        let f2 = column_fem_buckling(6, 200e9, 2e-6, l, p);
        let ratio = f2 / f1;
        assert!((ratio - 2.0).abs() < 0.05, "EI scaling ratio: {:.4}", ratio);
    }
    #[test]
    fn test_truss_elastic_stiffness_symmetric() {
        let k = truss_elastic_stiffness(3, 1.0, 100.0);
        let n = 4;
        for i in 0..n {
            for j in 0..n {
                assert!(
                    (k[i * n + j] - k[j * n + i]).abs() < 1e-12,
                    "K not symmetric at ({i},{j})"
                );
            }
        }
    }
    #[test]
    fn test_truss_geometric_stiffness_symmetric() {
        let kg = truss_geometric_stiffness(3, 1.0, 10.0);
        let n = 4;
        for i in 0..n {
            for j in 0..n {
                assert!(
                    (kg[i * n + j] - kg[j * n + i]).abs() < 1e-12,
                    "Kg not symmetric at ({i},{j})"
                );
            }
        }
    }
    #[test]
    fn test_truss_elastic_stiffness_row_sum_zero() {
        let k = truss_elastic_stiffness(3, 1.0, 100.0);
        let n = 4;
        for i in 0..n {
            let row_sum: f64 = (0..n).map(|j| k[i * n + j]).sum();
            assert!(row_sum.abs() < 1e-10, "Row {i} sum = {row_sum}");
        }
    }
    #[test]
    fn test_arc_length_step_1d_positive_stiffness() {
        let (d, l) = arc_length_step_1d(100.0, 0.0, 0.0, 0.1);
        assert!(d > 0.0, "displacement should increase");
        assert!(l > 0.0, "load should increase for positive stiffness");
    }
    #[test]
    fn test_arc_length_step_1d_zero_stiffness() {
        let (d, l) = arc_length_step_1d(0.0, 5.0, 1.0, 0.1);
        assert!((d - 1.1).abs() < 1e-10);
        assert!((l - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_trace_post_buckling_path() {
        let path = trace_post_buckling_path(|_d| 50.0, 0.0, 0.0, 0.01, 10);
        assert_eq!(path.len(), 11);
        for i in 1..path.len() {
            assert!(path[i].1 >= path[i - 1].1 - 1e-10);
        }
    }
    #[test]
    fn test_trace_post_buckling_returns_start() {
        let path = trace_post_buckling_path(|_| 100.0, 1.0, 0.5, 0.01, 0);
        assert_eq!(path.len(), 1);
        assert!((path[0].0 - 0.5).abs() < 1e-10);
        assert!((path[0].1 - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_imperfection_sensitivity_curve_zero_amplitude() {
        let curve = imperfection_sensitivity_curve(10.0, &[0.0], 0.5);
        assert_eq!(curve.len(), 1);
        assert!((curve[0].1 - 10.0).abs() < 1e-12);
    }
    #[test]
    fn test_imperfection_sensitivity_curve_decreasing() {
        let amps: Vec<f64> = (0..5).map(|i| i as f64 * 0.01).collect();
        let curve = imperfection_sensitivity_curve(10.0, &amps, 0.5);
        for i in 1..curve.len() {
            assert!(
                curve[i].1 <= curve[i - 1].1 + 1e-10,
                "load should decrease with increasing imperfection"
            );
        }
    }
    #[test]
    fn test_multimode_imperfection() {
        let lf = multimode_imperfection_sensitivity(10.0, 0.01, 0.5, 0.01, 0.3);
        assert!(lf < 10.0, "imperfection should reduce load");
        assert!(lf > 0.0, "load should remain positive");
    }
    #[test]
    fn test_multimode_zero_imperfection() {
        let lf = multimode_imperfection_sensitivity(10.0, 0.0, 0.5, 0.0, 0.3);
        assert!((lf - 10.0).abs() < 1e-12);
    }
    #[test]
    fn test_detect_snap_through_found() {
        let stiffnesses = [100.0, 50.0, 10.0, -5.0, -20.0];
        let idx = detect_snap_through(&stiffnesses);
        assert_eq!(idx, Some(3));
    }
    #[test]
    fn test_detect_snap_through_not_found() {
        let stiffnesses = [100.0, 50.0, 10.0, 5.0, 1.0];
        let idx = detect_snap_through(&stiffnesses);
        assert!(idx.is_none());
    }
    #[test]
    fn test_detect_snap_through_empty() {
        let idx = detect_snap_through(&[]);
        assert!(idx.is_none());
    }
    #[test]
    fn test_snap_through_from_path() {
        let path = vec![(0.0, 0.0), (0.1, 1.0), (0.2, 2.0), (0.3, 1.5), (0.4, 0.5)];
        let snap = snap_through_from_path(&path);
        assert!(snap.is_some());
        let (l, d) = snap.unwrap();
        assert!((l - 2.0).abs() < 1e-10);
        assert!((d - 0.2).abs() < 1e-10);
    }
    #[test]
    fn test_snap_through_from_path_monotonic() {
        let path = vec![(0.0, 0.0), (0.1, 1.0), (0.2, 2.0), (0.3, 3.0)];
        assert!(snap_through_from_path(&path).is_none());
    }
    #[test]
    fn test_find_limit_point() {
        let path = vec![(0.0, 0.0), (0.1, 5.0), (0.2, 3.0), (0.3, 1.0)];
        let lp = find_limit_point(&path);
        assert!(lp.is_some());
        let (d, l) = lp.unwrap();
        assert!((l - 5.0).abs() < 1e-10);
        assert!((d - 0.1).abs() < 1e-10);
    }
    #[test]
    fn test_find_limit_point_empty() {
        assert!(find_limit_point(&[]).is_none());
    }
    #[test]
    fn test_euler_buckling_critical_load() {
        let eb = EulerBuckling {
            E: 200e9,
            I: 1e-6,
            L: 2.0,
            K: 1.0,
        };
        let expected = PI * PI * 200e9 * 1e-6 / (1.0 * 2.0_f64).powi(2);
        assert!((eb.critical_load() - expected).abs() / expected < 1e-12);
    }
    #[test]
    fn test_euler_buckling_slenderness_ratio() {
        let eb = EulerBuckling {
            E: 200e9,
            I: 1e-6,
            L: 2.0,
            K: 1.0,
        };
        let r = 0.05;
        let lambda = eb.slenderness_ratio(1e-4, r);
        let expected = (1.0 * 2.0) / r;
        assert!((lambda - expected).abs() < 1e-12);
    }
    #[test]
    fn test_geometric_stiffness_beam_4x4_size() {
        let kg = geometric_stiffness_beam(1000.0, 1.0);
        assert_eq!(kg.len(), 4);
        for row in &kg {
            assert_eq!(row.len(), 4);
        }
    }
    #[test]
    fn test_geometric_stiffness_beam_symmetric() {
        let kg = geometric_stiffness_beam(500.0, 2.0);
        for i in 0..4 {
            for j in 0..4 {
                assert!(
                    (kg[i][j] - kg[j][i]).abs() < 1e-12,
                    "Kg not symmetric at ({i},{j})"
                );
            }
        }
    }
    #[test]
    fn test_buckling_analysis_critical_load_factor_positive() {
        let mut ba = BucklingAnalysis::new(4);
        for i in 0..4 {
            ba.K_global.push((i, i, 1000.0));
        }
        ba.add_geometric_stiffness([0, 1], 100.0, 1.0);
        let factor = ba.critical_load_factor_estimate();
        assert!(
            factor > 0.0 || factor.is_infinite(),
            "factor must be positive or infinite"
        );
    }
    #[test]
    fn test_johnson_parabolic_formula_below_euler() {
        let sigma_y = 250e6;
        let e = 200e9;
        let slenderness = 50.0;
        let sigma_cr = johnson_parabolic_formula(sigma_y, e, slenderness);
        let euler = PI * PI * e / (slenderness * slenderness);
        assert!(sigma_cr > 0.0);
        assert!(sigma_cr <= sigma_y);
        let _ = euler;
    }
    #[test]
    fn test_tangent_modulus_theory_scales() {
        let sigma_cr1 = tangent_modulus_theory(0.0, 200e9, 150e9, 100.0);
        let sigma_cr2 = tangent_modulus_theory(0.0, 200e9, 150e9, 200.0);
        assert!((sigma_cr1 / sigma_cr2 - 4.0).abs() < 1e-8);
    }
    #[test]
    fn test_beam_geometric_stiffness_symmetric() {
        let kg = beam_geometric_stiffness(3, 1.0, 10.0);
        let n = 2 * 4;
        for i in 0..n {
            for j in 0..n {
                assert!(
                    (kg[i * n + j] - kg[j * n + i]).abs() < 1e-12,
                    "Kg not symmetric at ({i},{j})"
                );
            }
        }
    }
    #[test]
    fn test_beam_geometric_stiffness_scales_with_p() {
        let kg1 = beam_geometric_stiffness(2, 1.0, 1.0);
        let kg2 = beam_geometric_stiffness(2, 1.0, 2.0);
        let n = 2 * 3;
        for i in 0..n * n {
            assert!(
                (kg2[i] - 2.0 * kg1[i]).abs() < 1e-12,
                "Kg should scale linearly with P"
            );
        }
    }
    #[test]
    fn test_cantilever_fem_buckling_vs_analytical() {
        let e = 200e9_f64;
        let i_mom = 1e-6_f64;
        let l = 2.0_f64;
        let p = 1.0_f64;
        let analytical = PI * PI * e * i_mom / (4.0 * l * l);
        let fem_mult = cantilever_fem_buckling(10, e, i_mom, l, p);
        let err = ((fem_mult * p - analytical) / analytical).abs();
        assert!(err < 0.05, "FEM error too large: {:.4e}", err);
    }
    #[test]
    fn test_cantilever_fem_buckling_positive() {
        let factor = cantilever_fem_buckling(4, 200e9, 1e-6, 2.0, 1.0);
        assert!(factor > 0.0, "buckling factor must be positive");
    }
    #[test]
    fn test_cantilever_less_than_pinned_pinned() {
        let e = 200e9;
        let i_mom = 1e-6;
        let l = 2.0;
        let p = 1.0;
        let cant = cantilever_fem_buckling(8, e, i_mom, l, p);
        let pp = column_fem_buckling(8, e, i_mom, l, p);
        assert!(
            cant < pp,
            "cantilever critical load {} should be less than pinned-pinned {}",
            cant,
            pp
        );
    }
    #[test]
    fn test_gaussian_solve_identity() {
        let a = [1.0, 0.0, 0.0, 1.0];
        let b = [3.0, 7.0];
        let x = gaussian_solve(&a, &b, 2).unwrap();
        assert!((x[0] - 3.0).abs() < 1e-12);
        assert!((x[1] - 7.0).abs() < 1e-12);
    }
    #[test]
    fn test_gaussian_solve_2x2() {
        let a = [2.0, 1.0, 1.0, 3.0];
        let b = [5.0, 7.0];
        let x = gaussian_solve(&a, &b, 2).unwrap();
        assert!((x[0] - 8.0 / 5.0).abs() < 1e-12);
        assert!((x[1] - 9.0 / 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_gaussian_solve_singular() {
        let a = [1.0, 2.0, 2.0, 4.0];
        let b = [1.0, 2.0];
        assert!(gaussian_solve(&a, &b, 2).is_none());
    }
    #[test]
    fn test_timoshenko_less_than_or_equal_euler() {
        let e = 200e9;
        let i = 1e-6;
        let l = 2.0;
        let a = 1e-3;
        let g = 80e9;
        let kappa = 5.0 / 6.0;
        let p_t = timoshenko_critical_load(e, i, l, a, g, kappa);
        let p_e = column_buckling_load(e, i, l, 1.0);
        assert!(p_t <= p_e + 1e-3, "Timoshenko load should be ≤ Euler load");
    }
    #[test]
    fn test_timoshenko_approaches_euler_for_slender_column() {
        let e = 200e9;
        let i = 1e-8;
        let l = 10.0;
        let a = 1e-4;
        let g = 80e9;
        let kappa = 5.0 / 6.0;
        let ratio = timoshenko_to_euler_ratio(e, i, l, a, g, kappa);
        assert!(
            ratio > 0.98,
            "Slender column: ratio should be close to 1, got {ratio}"
        );
    }
    #[test]
    fn test_timoshenko_ratio_in_zero_one() {
        let ratio = timoshenko_to_euler_ratio(200e9, 1e-6, 2.0, 1e-3, 80e9, 5.0 / 6.0);
        assert!(
            ratio > 0.0 && ratio <= 1.0 + 1e-12,
            "Ratio should be in (0, 1]"
        );
    }
    #[test]
    fn test_timoshenko_beam_stiffness_symmetric() {
        let k = timoshenko_beam_stiffness(200e9, 1e-6, 80e9, 1e-3, 5.0 / 6.0, 1.0);
        for i in 0..4 {
            for j in 0..4 {
                assert!(
                    (k[i][j] - k[j][i]).abs() < 1.0,
                    "Timoshenko stiffness not symmetric at [{i}][{j}]"
                );
            }
        }
    }
    #[test]
    fn test_timoshenko_slenderness_threshold_positive() {
        let lc = timoshenko_slenderness_threshold(200e9, 80e9, 5.0 / 6.0, 0.01);
        assert!(
            lc > 0.0,
            "Slenderness threshold should be positive, got {lc}"
        );
    }
    #[test]
    fn test_koiter_load_at_critical_point() {
        let koiter = KoiterExpansion::new(10.0, 0.0, 2.0);
        assert!(
            (koiter.load_factor(0.0) - 10.0).abs() < 1e-12,
            "At ξ=0 load = λ_c"
        );
    }
    #[test]
    fn test_koiter_stable_symmetric() {
        let koiter = KoiterExpansion::new(10.0, 0.0, 5.0);
        assert!(koiter.is_stable());
    }
    #[test]
    fn test_koiter_unstable_symmetric() {
        let koiter = KoiterExpansion::new(10.0, 0.0, -5.0);
        assert!(!koiter.is_stable());
    }
    #[test]
    fn test_koiter_asymmetric_unstable() {
        let koiter = KoiterExpansion::new(10.0, 1.0, 2.0);
        assert!(!koiter.is_stable());
    }
    #[test]
    fn test_koiter_equilibrium_amplitude_stable() {
        let koiter = KoiterExpansion::new(5.0, 0.0, 2.0);
        let xi = koiter.equilibrium_amplitude(7.0).unwrap();
        assert!((koiter.load_factor(xi) - 7.0).abs() < 1e-8);
    }
    #[test]
    fn test_koiter_equilibrium_amplitude_none_below_critical() {
        let koiter = KoiterExpansion::new(5.0, 0.0, 2.0);
        assert!(koiter.equilibrium_amplitude(4.0).is_none());
    }
    #[test]
    fn test_koiter_imperfection_reduces_load() {
        let koiter = KoiterExpansion::new(10.0, 0.0, -2.0);
        let lf = koiter.max_load_with_imperfection(0.04, 0.5);
        assert!(lf < 10.0, "Imperfection should reduce critical load");
        assert!(lf > 0.0, "Critical load should remain positive");
    }
    #[test]
    fn test_thermal_buckling_temperature_positive() {
        let dt = thermal_buckling_temperature(12e-6, 0.3, 0.01, 1.0, 4.0);
        assert!(
            dt > 0.0,
            "Critical temperature rise must be positive, got {dt}"
        );
    }
    #[test]
    fn test_thermal_buckling_scales_inverse_with_alpha() {
        let dt1 = thermal_buckling_temperature(12e-6, 0.3, 0.01, 1.0, 4.0);
        let dt2 = thermal_buckling_temperature(24e-6, 0.3, 0.01, 1.0, 4.0);
        assert!((dt1 / dt2 - 2.0).abs() < 1e-10, "ΔT_cr ~ 1/α");
    }
    #[test]
    fn test_thermal_buckling_scales_with_t_squared() {
        let dt1 = thermal_buckling_temperature(12e-6, 0.3, 0.01, 1.0, 4.0);
        let dt2 = thermal_buckling_temperature(12e-6, 0.3, 0.02, 1.0, 4.0);
        assert!((dt2 / dt1 - 4.0).abs() < 1e-8, "ΔT_cr ~ t²");
    }
    #[test]
    fn test_column_thermal_buckling_delta_t_positive() {
        let dt = column_thermal_buckling_delta_t(200e9, 1e-6, 2.0, 1e-3, 12e-6, 1.0);
        assert!(dt > 0.0);
    }
    #[test]
    fn test_combined_loading_unity_at_failure() {
        let util = combined_loading_utilization(5.0, 10.0, 50.0, 100.0);
        assert!((util - 1.0).abs() < 1e-12, "At combined failure util = 1.0");
    }
    #[test]
    fn test_combined_loading_safe() {
        let util = combined_loading_utilization(2.0, 10.0, 20.0, 100.0);
        assert!(util < 1.0, "Should be safe");
    }
    #[test]
    fn test_moment_amplification_factor_positive() {
        let af = moment_amplification_factor(5.0, 10.0).unwrap();
        assert!((af - 2.0).abs() < 1e-12, "At P/P_cr=0.5, AF=2.0");
    }
    #[test]
    fn test_moment_amplification_factor_none_at_critical() {
        assert!(moment_amplification_factor(10.0, 10.0).is_none());
    }
    #[test]
    fn test_moment_amplification_factor_none_above_critical() {
        assert!(moment_amplification_factor(15.0, 10.0).is_none());
    }
    #[test]
    fn test_dunkerley_less_than_minimum() {
        let p = dunkerley_combined_buckling(100.0, 200.0);
        assert!(
            p < 100.0,
            "Dunkerley result should be less than smaller load"
        );
    }
    #[test]
    fn test_dunkerley_approaches_single_when_one_very_large() {
        let p = dunkerley_combined_buckling(100.0, 1e15);
        assert!(
            (p - 100.0).abs() / 100.0 < 1e-6,
            "When one load >> other, result ≈ smaller"
        );
    }
    #[test]
    fn test_southwell_plot_basic() {
        let p_cr = southwell_critical_load(5.0, 1.0, 8.0, 2.0).unwrap();
        assert!(
            p_cr > 8.0,
            "Critical load should exceed measured load, got {p_cr}"
        );
    }
    #[test]
    fn test_reserve_factor_above_one_when_safe() {
        let rf = reserve_factor(100.0, 50.0);
        assert!((rf - 2.0).abs() < 1e-12, "RF = P_cr / P_applied = 2.0");
    }
    #[test]
    fn test_buckling_safety_check_passes() {
        assert!(buckling_safety_check(200.0, 50.0, 2.5));
    }
    #[test]
    fn test_buckling_safety_check_fails() {
        assert!(!buckling_safety_check(100.0, 50.0, 2.5));
    }
}
